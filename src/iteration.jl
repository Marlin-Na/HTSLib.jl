
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
