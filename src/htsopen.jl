
# Wrapper of htsFile with customized backend
# ------------------------------------------

## ? It seems that htsindex does not have to associate with htsFile?
## ? May consider removing sethtsindex! and use a new struct to manage index's lifetime.

export SamHeader
export HTSReadWriter

function remote_io(source::AbstractString; write=false, gcsclient = nothing)::IO
    if startswith(source, "gs://")
        if isnothing(gcsclient)
            gcsclient = HttpIO.GCSClient()
        end
        write == true && throw(ArgumentError("GCSFileIO does not support write"))
        return GCSFileIO(source, gcsclient)
    end
    if startswith(source, "http://") || startswith(source, "https://")
        write == true && throw(ArgumentError("HttpFileIO does not support write"))
        return HttpFileIO(source)
    end
    # assume local file
    mode = (write == true) ? "r+" : "r"
    return open(source, mode)
end

## SamHeader should be treated like an immutable object
mutable struct SamHeader
    sam_hdr::Ptr{htslib.sam_hdr_t}
    # below are cached fields
    # TODO: remove these cached fields?
    text::String
    lines::Vector{SubString}
    function SamHeader(ptr::Ptr)
        ptr == C_NULL && error("pointer is invalid")
        ans = new(ptr)
        Base.finalizer(ans) do ans
            ans.sam_hdr == C_NULL && error("pointer is invalid")
            htslib.sam_hdr_destroy(ans.sam_hdr)
        end
        ans.text = unsafe_string(htslib.sam_hdr_str(ans.sam_hdr))
        ans.lines = collect(filter(s -> s != "", split(ans.text, "\n")))
        ans
    end
end

function Base.pointer(x::SamHeader)
    x.sam_hdr
end

function SamHeader(x::SamHeader)
    dup = htslib.sam_hdr_dup(x.sam_hdr)
    dup == C_NULL && error("failed to duplicate sam_hdr: $(x.sam_hdr)")
    SamHeader(dup)
end

### TODO: wrap sam_hdr_nref, sam_hdr_name2tid, sam_hdr_tid2name and sam_hdr_tid2len

mutable struct HTSReadWriter{T<:IO}
    # hold the reference to objects below to avoid gc
    io::T
    hfile_backend::Ref{htslib.hFILE_backend}
    # htsFile pointer
    htsFile::Ptr{htslib.htsFile}
    # hFILE pointer
    hfile::Ptr{htslib.hFILE_julia}
    # samheader pointer
    samheader::SamHeader
    function HTSReadWriter{T}() where {T}
        new{T}()
    end
end

function Base.pointer(x::HTSReadWriter)
    x.htsFile
end

function BioGenerics.header(hf::HTSReadWriter)
    hf.samheader
end

function BioGenerics.IO.stream(hf::HTSReadWriter)
    hf.io
end

function SamHeader(hf::HTSReadWriter)
    SamHeader(header(hf)) # copy
end

function SamHeader(x::String)
    ptr = htslib.sam_hdr_parse(sizeof(x), pointer(x))
    ptr == C_NULL && error("error parsing samheader")
    SamHeader(ptr)
end

function Base.string(x::SamHeader)
    x.text
end

function Base.show(io::IO, x::SamHeader)
    s = string(x)
    lines = split(s, "\n")
    lines = filter(l -> l != "", lines)
    SQlines = filter(l -> l[1:3] == "@SQ", lines)
    first_seq = split(split(SQlines[1], "\t")[2], ":")[2]
    last_seq  = split(split(SQlines[end], "\t")[2], ":")[2]
    print(io, typeof(x), ": ", length(lines), " lines; ", length(SQlines), " sequences '$(first_seq)...$(last_seq)';")
    x
end

function HTSReadWriter(source::AbstractString; index=nothing, write=nothing, gcsclient=nothing, kwargs...)
    # reuse gcsclient to download bam index and construct the reader
    if startswith(source, "gs://") && isnothing(gcsclient)
        gcsclient = GCSClient()
    end
    if isa(index, AbstractString) && startswith(index, "gs://")
        index = read(remote_io(index; gcsclient=gcsclient))::Vector{UInt8}
    end

    bamio = remote_io(source; write=write, gcsclient=gcsclient)
    HTSReadWriter(bamio; index=index, write=write, kwargs...)
end

function HTSReadWriter(io::IO; index=nothing, write=nothing)
    write2 = false
    if isnothing(write) && applicable(iswritable, io) && iswritable(io)
        write2 = true
    end
    mode = write2 ? "r+" : "r"

    hfileptr = htslib.hfile_init(sizeof(htslib.hFILE_julia), pointer(mode), 0)
    hfileptr = convert(Ptr{htslib.hFILE_julia}, hfileptr)
    hfile_backend = Ref(htslib.hfile_julia_backend())
    ## Add argument's pointer to the struct and assign julia's backend implementation
    let
        hfile = unsafe_load(hfileptr)
        hfile_base2 = htslib.hFILE_base(
            hfile.base.buffer   ,
            hfile.base.beginz   ,
            hfile.base.endz     ,
            hfile.base.limit    ,
            pointer_from_objref(hfile_backend),
            hfile.base.offset   ,
            hfile.base.remaining,
        )
        hfile2 = htslib.hFILE_julia(hfile_base2, pointer_from_objref(io))
        unsafe_store!(hfileptr, hfile2)
    end
    ## htsFile pointer
    htsFileptr = htslib.hts_hopen(hfileptr, pointer("/dummy"), pointer(mode))
    ## wrapper type holding both julia object and pointers managed in C
    #h = HTSReadWriter{typeof(io)}(htsFileptr, hfileptr, io, hfile_backend)
    h = HTSReadWriter{typeof(io)}()
    let
        h.io = io
        h.hfile_backend = hfile_backend
        h.htsFile = htsFileptr
        h.hfile = hfileptr

        Base.finalizer(h) do h
            #h.hfile == C_NULL && error("invalid pointer")
            #h.htsFile == C_NULL && error("invalid pointer")
            if h.htsFile != C_NULL
                ret = htslib.hts_close(h.htsFile)
                ret < 0 && Base.systemerror("error running hts_close with exit code $ret")
            end
            # Looks like that hts_close will also free the underlying hfile. So we don't need it here.
            # htslib.hfile_destroy(hfileptr)
        end
    end
    # TODO
    # let
    #     # set number of threads to use
    #     hts_set_threads!(h, nthreads) 
    # end
    let
        # read sam header
        sam_hdr = htslib.sam_hdr_read(htsFileptr)
        sam_hdr == C_NULL && error("error reading sam header")
        h.samheader = SamHeader(sam_hdr)
    end
    # load index.
    htssetindex!(h, index)
    h
end

# TODO: also implement it for HTSRegionsIterator

function Base.open(::Type{HTSReadWriter}, args...; kwargs...)
    HTSReadWriter(args...; kwargs...)
end

function Base.open(f::Function, ::Type{HTSReadWriter}, args...; kwargs...)
    io = HTSReadWriter(HTSReadWriter, args...; kwargs...)
    try
        f(io)
    finally
        close(io)
    end
end

"""
    htssetindex!(htsio::HTSReadWriter, index::AbstractString)

Load index of hts file. Supports local file, http url and gs url.
"""
#function htssetindex!(htsio::HTSReadWriter, index)
#    htssetindex!(htsio, HTSIndexData(index))
#end

function htssetindex!(htsio::HTSReadWriter, index::AbstractString)
    # short circuit for local file (avoid current roundtrip of disk->memory->disk->memory)
    if index isa AbstractString && (~occursin(r"^.*://.*$", index))
        ptr = htslib.hts_idx_load2(pointer("/dummy"), pointer(index))
        ptr == C_NULL && error("error loading index file $(index)")
        let
            tmp = unsafe_load(pointer(htsio))
            tmp = @set tmp.idx = ptr
            unsafe_store!(pointer(htsio), tmp)
        end
        return htsio
    end

    htssetindex!(htsio, remote_io(index))
end

# remove index
function htssetindex!(htsio::HTSReadWriter, index::Nothing)
    # hts_close will also free the index, otherwise we need to call hts_idx_destroy
    oldidx = unsafe_load(pointer(htsio)).idx
    let
        tmp = unsafe_load(pointer(htsio))
        tmp = @set tmp.idx = C_NULL
        unsafe_store!(pointer(htsio), tmp)
    end
    oldidx != C_NULL && htslib.hts_idx_destroy(oldidx)
    return htsio
end

function htssetindex!(htsio::HTSReadWriter, index::IO)
    data = read(index)::Vector{UInt8}
    htssetindex!(htsio, data)
end

function htssetindex!(htsio::HTSReadWriter, index::Vector{UInt8})
    # this is a dirty solution
    ptr = mktemp() do tmppath, tmpio
        write(tmpio, index)
        flush(tmpio)
        ptr = GC.@preserve tmppath htslib.hts_idx_load2(pointer("/dummy"), pointer(tmppath))
        ptr == C_NULL && throw("Failed to load hts index file")
        ptr
    end

    # if index already exists, call hts_idx_destroy on the old index
    oldidx = unsafe_load(pointer(htsio)).idx

    # assign new index
    let
        tmp = unsafe_load(pointer(htsio))
        tmp = @set tmp.idx = ptr
        unsafe_store!(pointer(htsio), tmp)
    end

    oldidx != C_NULL && htslib.hts_idx_destroy(oldidx)
    return htsio
end

function Base.show(io::IO, hf::HTSReadWriter)
    if pointer(hf) == C_NULL
        print(io, typeof(hf), " Invalid HTSReadWriter")
        return nothing
    end
    nseq = length(filter(startswith("@SQ"), split(string(header(hf)), "\n")))
    print(io, typeof(hf), " ", unsafe_string(htslib.hts_format_description(htslib.hts_get_format(hf.htsFile))), " with ", nseq, " reference sequences")
end

# Methods
# -------

# During initialization

function hts_set_threads!(hf::HTSReadWriter, n)
    ret = htslib.hts_set_threads(hf.htsFile, convert(Int32, n))
    ret != 0 && error("error running hts_set_threads")
    hf
end

function hts_set_cache_size!(hf::HTSReadWriter, n)
    htslib.hts_set_cache_size(hf.htsFile, convert(Int32, n))
end

# Read

function BioGenerics.IO.tryread!(hf::HTSReadWriter, record::HTSRecord)::Union{HTSRecord,Nothing}
    ret = htslib.sam_read1(pointer(hf), pointer(hf.samheader), record.ptr)
    ret < -1 && error("error reading record")
    ret == -1 && return nothing
    return record
end

function Base.read!(hf::HTSReadWriter, record::HTSRecord)
    isnothing(tryread!(hf, record)) && throw(EOFError())
    record
end

function Base.read(hf::HTSReadWriter)
    read!(hf, HTSRecord())
end

# eof, flush, close

# TODO: how to implement eof and flush?

function Base.close(hf::HTSReadWriter)
    # finalizer will run hts_close, which will close the underlying IO from julia
    Base.finalize(hf)
    hf.htsFile = C_NULL
    return nothing
end

# write

# TODO
