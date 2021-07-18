
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

