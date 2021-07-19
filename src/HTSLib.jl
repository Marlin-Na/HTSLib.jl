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
    header,
    IO.tryread!,
    IO.stream

# BioGenerics reexport
export seqname,
    sequence,
    metadata,
    leftposition,
    rightposition,
    isoverlapping,
    metainfoval,
    metainfotag,
    header,
    tryread!,
    stream

export sequence!,
    setsequence!,
    seqlength,
    mappingquality,
    refid,
    refname,
    quality,
    quality!,
    setquality!,
    queryname

export htslib
export eachrecord

using Printf: @printf
using Setfield: @set

import HttpIO: HttpIO, PoorGCloudAuth, GCSClient, GCSFileIO, HttpFileIO

import TranscodingStreams

include("htslib/htslib.jl")
include("bamrecord.jl")
include("htsopen.jl")
include("iteration.jl")

"""
The version of htslib.
"""
const HTSLIB_VERSION = VersionNumber(unsafe_string(htslib.hts_version()))


end # module
