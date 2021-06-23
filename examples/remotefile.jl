
using HTSLib
import Downloads

"""
    count_remote_file_records(url::AbstractString)

Iterate over the bam file. Use `curl` or `gsutil` to stream the data.
"""
function count_remote_file_records(url::AbstractString)
    if startswith(url, "http://") || startswith(url, "https://")
        stream = open(pipeline(`curl --silent $url`), "r")
    elseif startswith(url, "gs://")
        stream = open(pipeline(`gsutil cat $url`), "r")
    else
        # assume this is a local file
        stream = open(url, "r")
    end

    ## htsopen may use arbitrary IO object in Julia
    count_remote_file_records(stream)
end

function count_remote_file_records(io::IO)
    nreads::Int = 0
    for record in HTSLib.htsopen(io, "r")
        nreads += 1
    end
    println("number of reads: $nreads")
end

# sample test files:
#   - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam
#   - gs://broad-public-datasets/NA12878/NA12878.cram

if PROGRAM_FILE == relpath(@__FILE__)
    @assert length(ARGS) == 1
    url = ARGS[1]
    count_remote_file_records(url)
end
