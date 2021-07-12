
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

function BamRecord(x::BamRecord)
    copy(x)
end

function Base.copy!(dst::BamRecord, src::BamRecord)
    GC.@preserve src begin
        p = htslib.bam_copy1(pointer(dst), pointer(src))
        p == C_NULL && error("error copying bam record $src")
    end
    dst
end

function Base.copy(record::BamRecord)
    ans = BamRecord()
    copy!(ans, record)
    ans
end

function Base.pointer(x::BamRecord)
    x.ptr
end

function Base.show(io::IO, record::BamRecord)
    GC.@preserve record begin
        ptr = record.ptr
        ptr == C_NULL && error("pointer is invalid")
        bam = unsafe_load(ptr)
        @printf(io, "%s(<%p>) %d:%d-%d ", summary(record), ptr, seqlevel(record), leftposition(record), rightposition(record))
    end
    print(io, "\n")
    print(io, "queryname: $(queryname(record))\n")
    print(io, "mapqual  : $(mappingquality(record))\n")
    print(io, "sequence : ")
    for s in sequence(record)
        print(s)
    end
    print(io, "\n")
    print(io, "quality  : ")
    for q in quality(record)
        print(io, quality_char(q))
    end
    record
end

function quality_char(qual::UInt8)
    # taken from jakobnissen
	if qual < 10
		return  ' '
	end
	if qual < 15
		return '▁'
	end
	if qual < 20
		return '▂'
	end
	if qual < 25
		return '▃'
	end
	if qual < 30
		return '▄'
	end
	if qual < 35
		return '▆'
	end
	if qual < 40
		return '▇'
	end
	if qual < 255
		return '█'
	end
	return '?'
end

@inline function Base.getproperty(view::BamRecord, name::Symbol)
    GC.@preserve view begin
        bam = unsafe_load(getfield(view, :ptr))
    end
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


##### Accessor functions

const seq_nt16_str = ('=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N')

function BioGenerics.leftposition(record::BamRecord)::Int64
    record.pos + 1 # start pos is zero-based
end

function BioGenerics.rightposition(record::BamRecord)::Int64
    GC.@preserve record begin
        htslib.bam_endpos(pointer(record))
    end
end

function seqlength(record::BamRecord)::Int32
    record.l_qseq
end

function BioGenerics.sequence(record::BamRecord, i::Integer)::Char
    GC.@preserve record begin
        seq_ptr = htslib.bam_get_seq(pointer(record))
        idx = htslib.bam_seqi(seq_ptr, i-UInt8(1)) + 1
        seq_nt16_str[idx]
    end
end

function sequence!(record::BamRecord, array::AbstractVector{Char})
    seql = seqlength(record)
    empty!(array)
    sizehint!(array, seql)
    GC.@preserve record begin
        seq_ptr = htslib.bam_get_seq(pointer(record))
        for i in 1:seql
            push!(array, seq_nt16_str[htslib.bam_seqi(seq_ptr, i-UInt8(1)) + 1])
        end
    end
    array
end

function BioGenerics.sequence(record::BamRecord)
    ans = Vector{Char}(undef, seqlength(record))
    sequence!(record, ans)
end

function mappingquality(record::BamRecord)::UInt8
    record.qual
end

function seqlevel(record::BamRecord)::Int32
    record.tid + Int32(1)
end

function quality!(record::BamRecord, array::AbstractVector)
    seql = seqlength(record)
    empty!(array)
    sizehint!(array, seql)
    GC.@preserve record begin
        qual_ptr = htslib.bam_get_qual(pointer(record))
        for i in 1:seql
            push!(array, unsafe_load(qual_ptr, i))
        end
    end
    array
end

function quality(record::BamRecord)
    #### TODO: quality score 255 represents missing. Do we want to return missing?
    ans = Vector{UInt8}(undef, seqlength(record))
    quality!(record, ans)
end

function setquality!(record::BamRecord, i::Integer, value::Integer)
    #### TODO: quality score 255 represents missing.?
    value = convert(UInt8, value)
    GC.@preserve record begin
        ## bounds check
        if i <= 0 || i > seqlength(record)
            throw(BoundsError(record, i))
        end
        qual_ptr = htslib.bam_get_qual(pointer(record))
        unsafe_store!(qual_ptr, value, i)
    end
    record
end

function setquality!(record::BamRecord, value::AbstractVector{<:Integer})
    seqlength(record) != length(value) && throw(BoundsError(value))
    for (i, v) in enumerate(value)
        setquality!(record, i, v)
    end
    record
end

function setsequence!(record::BamRecord, i::Integer, value::Char)
    GC.@preserve record begin
        ## bounds check
        if i <= 0 || i > seqlength(record) 
            throw(BoundsError(record, i))
        end
        idx = findfirst(x -> x == value, seq_nt16_str)
        isnothing(idx) && throw(DomainError(value))
        base_code = UInt8(idx - 1)
        htslib.bam_set_seqi(htslib.bam_get_seq(pointer(record)), i-1, base_code)
    end
    record
end

function queryname(record::BamRecord)::String
    # Note: l_qname is length of the name _plus_ number of NULs (1-4) for memory alignment
    GC.@preserve record begin
        l_qname = record.l_qname
        l_qname <= 0 && return "*" # If set to 0, the placeholder query name "*" will be used.
        cstr = htslib.bam_get_qname(pointer(record))
        unsafe_string(pointer(cstr))
    end
end
