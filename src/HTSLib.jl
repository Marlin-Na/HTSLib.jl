module HTSLib

export htslib
#export eachrecord

using Printf: @printf
using Setfield: @set

import TranscodingStreams

include("htslib/htslib.jl")

"""
The version of htslib.
"""
const HTSLIB_VERSION = VersionNumber(unsafe_string(htslib.hts_version()))


# Record view
# -----------

export BamRecord

"""
A light-weight wrapper of records.
"""
mutable struct BamRecord
    ptr::Ptr{htslib.bam1_t}
    function BamRecord()
        ptr = htslib.bam_init1()
        @assert ptr != C_NULL
        record = new(ptr)
        Base.finalizer(record) do record
            ptr = getfield(record, :ptr)
            ptr == C_NULL && error("pointer is invalid")
            htslib.bam_destroy1(ptr)
        end
        record
    end
end

function Base.pointer(x::BamRecord)
    x.ptr
end

function Base.show(io::IO, record::BamRecord)
    ptr = record.ptr
    ptr == C_NULL && error("pointer is invalid")
    bam = unsafe_load(ptr)
    @printf(io, "%s(<%d:%d@%p>) ", summary(record), bam.core.tid, bam.core.pos, ptr)
    record
end

@inline function Base.getproperty(view::BamRecord, name::Symbol)
    bam = unsafe_load(getfield(view, :ptr))
    if name == :pos
        return bam.core.pos
    elseif name == :tid
        return bam.core.tid
    elseif name == :bin
        return bam.core.bin
    elseif name == :qual
        return bam.core.qual
    elseif name == :l_extranul
        return bam.core.l_extranul
    elseif name == :flag
        return bam.core.flag
    elseif name == :l_qname
        return bam.core.l_qname
    elseif name == :n_cigar
        return bam.core.n_cigar
    elseif name == :l_qseq
        return bam.core.l_qseq
    elseif name == :mtid
        return bam.core.mtid
    elseif name == :mpos
        return bam.core.mpos
    elseif name == :isize
        return bam.core.isize
    else
        return getfield(view, name)
    end
end

include("htsopen.jl")


end # module
