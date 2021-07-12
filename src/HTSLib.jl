module HTSLib

import BioGenerics:
    BioGenerics,
    seqname,
    sequence,
    metadata,
    leftposition,
    rightposition,
    isoverlapping,
    metainfoval,
    metainfotag,
    header

# BioGenerics reexport
export seqname,
    sequence,
    metadata,
    leftposition,
    rightposition,
    isoverlapping,
    metainfoval,
    metainfotag,
    header

export sequence!,
    seqlength,
    mappingquality,
    seqlevel,
    quality,
    setquality!,
    setsequence!

export htslib
#export eachrecord

using Printf: @printf
using Setfield: @set

import TranscodingStreams

include("htslib/htslib.jl")
include("bamrecord.jl")
include("htsopen.jl")

"""
The version of htslib.
"""
const HTSLIB_VERSION = VersionNumber(unsafe_string(htslib.hts_version()))


end # module
