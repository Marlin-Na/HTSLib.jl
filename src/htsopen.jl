
# Wrapper of htsFile with customized backend
# ------------------------------------------

export SamHeader
export HTSIOWrapper

mutable struct SamHeader
    sam_hdr::Ptr{htslib.sam_hdr_t}
    # below are cached fields
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

mutable struct HTSIOWrapper{T<:IO}
    # hold the reference to objects below to avoid gc
    io::T
    hfile_backend::Ref{htslib.hFILE_backend}
    # htsFile pointer
    htsFile::Ptr{htslib.htsFile}
    # hFILE pointer
    hfile::Ptr{htslib.hFILE_julia}
    # samheader pointer
    samheader::SamHeader
    function HTSIOWrapper{T}() where {T}
        new{T}()
    end
end

function Base.pointer(x::HTSIOWrapper)
    x.htsFile
end

function SamHeader(hf::HTSIOWrapper)
    SamHeader(hf.samheader)
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

function htsopen(source::AbstractString; kwargs...)
    htsopen(remote_io(source); kwargs...)
end

function htsopen(io::IO; index=nothing, write=false)
    mode = write ? "r+" : "r"
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
    #h = HTSIOWrapper{typeof(io)}(htsFileptr, hfileptr, io, hfile_backend)
    h = HTSIOWrapper{typeof(io)}()
    let
        h.io = io
        h.hfile_backend = hfile_backend
        h.htsFile = htsFileptr
        h.hfile = hfileptr

        Base.finalizer(h) do h
            h.htsFile == C_NULL && error("invalid pointer")
            h.hfile == C_NULL && error("invalid pointer")
            ret = htslib.hts_close(h.htsFile)
            ret < 0 && Base.systemerror("error running hts_close with exit code $ret")
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

"""
    htssetindex!(htsio::HTSIOWrapper, index::AbstractString)

Load index of hts file. Supports local file, http url and gs url.
"""
function htssetindex!(htsio::HTSIOWrapper, index)
    htssetindex!(htsio, HTSIndexData(index))
end

function htssetindex!(htsio::HTSIOWrapper, index::AbstractString)
    # short circuit for local file (avoid current roundtrip of disk->memory->disk->memory)
    if index isa AbstractString && (~occursin(r"^.*://.*$", index))
        ptr = htslib.hts_idx_load2(pointer("/dummy"), pointer(index))
        ptr == C_NULL && error("error loading index file $(index)")
        let
            tmp = unsafe_load(pointer(htsio))
            @set tmp.idx = ptr
            unsafe_store!(pointer(htsio), tmp)
        end
        return htsio
    end

    htssetindex!(htsio, HTSIndexData(index))
end

function htssetindex!(htsio::HTSIOWrapper, index::Nothing)
    # hts_close will also free the index, otherwise we need to call hts_idx_destroy
    oldidx = unsafe_load(pointer(htsio)).idx
    let
        tmp = unsafe_load(pointer(htsio))
        @set tmp.idx = C_NULL
        unsafe_store!(pointer(htsio), tmp)
    end
    oldidx != C_NULL && htslib.hts_idx_destroy(oldidx)
    return htsio
end

function htssetindex!(htsio::HTSIOWrapper, index::HTSIndexData)
    content = parent(index)::Vector{UInt8}
    # this is a dirty solution
    ptr = mktemp() do tmppath, tmpio
        write(tmpio, content)
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
        @set tmp.idx = ptr
        unsafe_store!(pointer(htsio), tmp)
    end

    oldidx != C_NULL && htslib.hts_idx_destroy(oldidx)
    return htsio
end

function Base.show(io::IO, hf::HTSIOWrapper)
    nseq = length(filter(startswith("@SQ"), split(string(SamHeader(hf)), "\n")))
    print(io, typeof(hf), " ", unsafe_string(htslib.hts_format_description(htslib.hts_get_format(hf.htsFile))), " with ", nseq, " reference sequences")
end

# Methods
# -------

# During initialization

function hts_set_threads!(hf::HTSIOWrapper, n)
    ret = htslib.hts_set_threads(hf.htsFile, convert(Int32, n))
    ret != 0 && error("error running hts_set_threads")
    hf
end

function hts_set_cache_size!(hf::HTSIOWrapper, n)
    htslib.hts_set_cache_size(hf.htsFile, convert(Int32, n))
end

function Base.read!(hf::HTSIOWrapper, record::BamRecord)
    ret = htslib.sam_read1(hf.htsFile, hf.samheader.sam_hdr, record.ptr)
    ret < -1 && error("error reading record")
    ret == -1 && throw(EOFError())
    record
end

# Iterator
# --------

function Base.iterate(hf::HTSIOWrapper, state::Missing)
    record = BamRecord()
    try
        read!(hf, record)
    catch err
        if err isa EOFError
            return nothing
        else
            rethrow()
        end
    end
    return (record, state)
end

function Base.iterate(hf::HTSIOWrapper)
    state = missing
    iterate(hf, state)
end
