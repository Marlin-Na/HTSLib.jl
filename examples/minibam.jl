
using HTSLib
using HttpIO

"""
Generate minibam files for the given regions.
"""
function make_minibam(bamfile, indexfile, regions, outfile; gcsclient=nothing)
    if isnothing(gcsclient)
        gcsclient = GCSClient()
    end

    ## First pass to find all mate regions
    reader = HTSReadWriter(bamfile; index=indexfile, gcsclient=gcsclient)
    mate_regions = String[]
    for r in HTSRegionsIterator(reader, regions)
        reg = find_mate_region(header(reader), r)
        (!isnothing(reg)) && push!(mate_regions, reg)
    end

    ## Second pass to fetch reads and their mates
    all_regions = vcat(regions, mate_regions)
    writer = HTSReadWriter(open(outfile, "w"), "wbz")
    write(writer, header(reader)) # write header
    for r in HTSRegionsIterator(reader, all_regions)
        compress_qual!(r)
        write(writer, r) # write record
    end
    close(writer)
end

function find_mate_region(hdr, record; window=30)::Union{Nothing,String}
    if record.flag.paired
        ntid = next_refid(record)
        npos = next_refpos(record)
        if (!ismissing(ntid)) && (!ismissing(npos))
            nchr = unsafe_string(htslib.sam_hdr_tid2name(pointer(hdr), Cint(ntid - 1)))
            return "$nchr:$(Int(npos-window))-$(Int(npos+window))"
        end
    end
    return nothing
end

"""
Compress a record's quality score according to illumina scheme:
https://www.illumina.com/documents/products/whitepapers/whitepaper_datacompression.pdf
"""
function compress_qual!(record::HTSRecord)
    qual = quality(record)
    for i in 1:length(qual)
        new_score = ILLUMINA_8_BINNING_DICT[qual[i]]
        setquality!(record, i, new_score)
    end
end

const ILLUMINA_8_BINNING_DICT = Dict{UInt8, UInt8}()
ILLUMINA_8_BINNING_DICT[0] = 0
ILLUMINA_8_BINNING_DICT[1] = 1
for i in 2:9
    ILLUMINA_8_BINNING_DICT[i] = 6
end
for i in 10:19
    ILLUMINA_8_BINNING_DICT[i] = 15
end
for i in 20:24
    ILLUMINA_8_BINNING_DICT[i] = 22
end
for i in 25:29
    ILLUMINA_8_BINNING_DICT[i] = 27
end
for i in 30:34
    ILLUMINA_8_BINNING_DICT[i] = 33
end
for i in 35:39
    ILLUMINA_8_BINNING_DICT[i] = 37
end
for i in 40:254
    ILLUMINA_8_BINNING_DICT[i] = 40
end
ILLUMINA_8_BINNING_DICT[0xff] = 0xff


#### test 

# bam = "gs://fc-secure-66f5eeb9-27c4-4e5c-b9d6-0519aca5889d/7d418b8a-5b21-443e-8d9b-93081096dc7a/gdc_api_file_download/7423815f-7175-4716-ad1f-af4635994334/call-download_file/d5011161-3d6b-4802-8a64-614f6d2298f9_wgs_gdc_realn.bam"
# bai = "gs://fc-secure-66f5eeb9-27c4-4e5c-b9d6-0519aca5889d/3f2e4ec2-5eac-4ee2-9180-6896f354b423/gdc_api_file_download/f9b78731-6131-4d92-ba5d-29ed7dd5cabf/call-download_file/d5011161-3d6b-4802-8a64-614f6d2298f9_wgs_gdc_realn.bai"
# 
# make_minibam(bam, bai, ["chrX:101254500-101254600"], "/tmp/test.out.mini.bam")
