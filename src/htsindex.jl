# Read HTSIndex from different sources.

struct HTSIndexData
    data::Vector{UInt8}
end

function HTSIndexData(source::AbstractString)
    HTSIndexData(read_remote_data(source))
end

function HTSIndexData(source::IO)
    data = read(source)::Vector{UInt8}
    HTSIndexData(data)
end

function Base.parent(idxdata::HTSIndexData)
    idxdata.data
end

"""
Read data from a url or local file path.
"""
function read_remote_data(source::AbstractString)::Vector{UInt8}
    io = remote_io(source)
    return read(io)
end

# this is type unstable
function remote_io(source::AbstractString)::IO
    if startswith(source, "gs://")
        return GCSFileIO(source)
    end
    if startswith(source, "http://") || startswith(source, "https://")
        return HttpFileIO(index)
    end
    # assume local file
    return open(source, "r+")
end
