
# Iterator
# --------

function Base.eltype(::HTSReadWriter)
    HTSRecord
end

function Base.IteratorSize(::Type{T}) where T<: HTSReadWriter
    Base.SizeUnknown()
end

function Base.iterate(hf::HTSReadWriter, state::Missing)
    record = HTSRecord()
    isnothing(tryread!(hf, record)) && return nothing
    return (record, state)
end

function Base.iterate(hf::HTSReadWriter)
    state = missing
    iterate(hf, state)
end

#function eachrecord(hf::HTSReadWriter)
#end



# MultiRegion Iterator
# --------------------

export HTSRegionsIterator

"""
Iterate alignment records over the specified regions. If a read
overlaps more than one region, it will only be returned once.

The region string can take one of the following forms:

| region          | Outputs                                                     |
|:--------------- |:----------------------------------------------------------- |
| REF             | All reads with RNAME REF                                    |
| REF:            | All reads with RNAME REF                                    |
| REF:START       | Reads with RNAME REF overlapping START to end of REF        |
| REF:-END        | Reads with RNAME REF overlapping start of REF to END        |
| REF:START-END   | Reads with RNAME REF overlapping START to END               |
| .               | All reads from the start of the file                        |
| *               | Unmapped reads at the end of the file (RNAME '*' in SAM)    |

"""
HTSRegionsIterator

mutable struct HTSRegionsIterator{T}
    htsreader::HTSReadWriter{T}
    itr_ptr::Ptr{htslib.hts_itr_t}
    function HTSRegionsIterator(hf::HTSReadWriter{T}, regions::AbstractVector) where T
        if !isa(regions, Vector{String})
            return HTSRegionsIterator(hf, map(r -> region_string(hf, r), regions))
        end
        # also need to preserve pointer array from gc
        @assert isa(regions, Vector{String})
        regarray = pointer.(regions)
        GC.@preserve regions regarray begin
            idxptr = unsafe_load(pointer(hf)).idx
            hdrptr = pointer(header(hf))
            regarrayptr = convert(Ptr{Ptr{Char}}, pointer(regarray))
            itrptr = htslib.sam_itr_regarray(idxptr, hdrptr, regarrayptr, Cuint(length(regions)))
            ans = new{T}(hf, itrptr)
            Base.finalizer(ans) do ans
                ptr = pointer(ans)
                ptr != C_NULL && htslib.sam_itr_destroy(ptr)
                return nothing
            end
            ans
        end
    end
end

function HTSRegionsIterator(hf::HTSReadWriter, chr, start, stop)
    regions = String[]
    if isa(chr, AbstractString)
        chr = Ref(chr)
    end
    for (c, s, e) in zip(chr, start, stop)
        region = region_string(hf, Tuple{typeof(c),typeof(s),typeof(e)}((c, s, e)))
        push!(regions, region)
    end
    HTSRegionsIterator(hf, regions)
end

function HTSRegionsIterator(hf::HTSReadWriter, region::Union{Tuple, AbstractString})
    regions = Vector{typeof(region)}(undef, 1)
    regions[1] = region
    return HTSRegionsIterator(hf, regions)
end

# bad solution, this should be a feature request to htslib
function region_string(hf::HTSReadWriter, regionspec::Tuple{<:AbstractString,<:Real,<:Real})::String
    (chr, start, stop) = regionspec
    @assert start >= 0
    @assert stop >= 0
    start = Int(start)
    stop = Int(stop)
    return "$chr:$start-$stop"
end

function region_string(hf::HTSReadWriter, regionspec::Tuple{<:Real,<:Real,<:Real})::String
    (chr, start, stop) = regionspec
    @assert chr >= 1
    @assert start >= 0
    @assert stop >= 0
    # TODO: provide wrapper for sam_hdr_tid2name
    chr_str = unsafe_string(htslib.sam_hdr_tid2name(pointer(header(hf)), Cint(chr-1)))
    start = Int(start)
    stop = Int(stop)
    return "$chr_str:$start-$stop"
end

function region_string(hf::HTSReadWriter, regionspec::AbstractString)::String
    String(regionspec)
end

function BioGenerics.IO.tryread!(itr::HTSRegionsIterator, record::HTSRecord)::Union{HTSRecord,Nothing}
    GC.@preserve itr begin
        res = htslib.sam_itr_next(pointer(parent(itr)), pointer(itr), pointer(record))
        res < -1 && error("error reading record")
        res == -1 && return nothing
        return record
    end
end

function Base.pointer(x::HTSRegionsIterator)
    x.itr_ptr
end

function Base.parent(x::HTSRegionsIterator)
    x.htsreader
end

function Base.eltype(::HTSRegionsIterator)
    HTSRecord
end

function Base.IteratorSize(::Type{T}) where T <: HTSRegionsIterator
    Base.SizeUnknown()
end

function Base.iterate(itr::HTSRegionsIterator, state::Missing)
    res = tryread!(itr, HTSRecord())
    isnothing(res) && return nothing
    return (res, missing)
end

function Base.iterate(itr::HTSRegionsIterator)
    iterate(itr, missing)
end

## implementation of read! and read are common between HTSRegionsIterator and HTSReadWriter
function Base.read!(hf::HTSRegionsIterator, record::HTSRecord)
    isnothing(tryread!(hf, record)) && throw(EOFError())
    record
end

function Base.read(hf::HTSRegionsIterator)
    read!(hf, HTSRecord())
end
