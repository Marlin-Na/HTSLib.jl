
const HTSRecordFlagBits = (
    # the read is paired in sequencing, no matter whether it is mapped in a pair
    paired=UInt16(1),
    # the read is mapped in a proper pair
    proper_pair=UInt16(2),
    # the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
    unmap=UInt16(4),
    # the mate is unmapped
    munmap=UInt16(8),
    # the read is mapped to the reverse strand
    reverse=UInt16(16),
    # the mate is mapped to the reverse strand
    mreverse=UInt16(32),
    # this is read1
    read1=UInt16(64),
    # this is read2
    read2=UInt16(128),
    # not primary alignment
    secondary=UInt16(256),
    # QC failure
    qcfail=UInt16(512),
    # optical or PCR duplicate
    dup=UInt16(1024),
    # supplementary alignment
    supplementary=UInt16(2048),
)

"""
Access HTSRecord bitwise flags.
"""
struct HTSRecordFlag
    value::UInt16
end

function Base.convert(::Type{UInt16}, flag::HTSRecordFlag)
    getfield(flag, :value)
end

@inline function Base.getproperty(flag::HTSRecordFlag, name::Symbol)
    bit = getproperty(HTSRecordFlagBits, name)
    (convert(UInt16, flag) & bit) > 0
end

function Base.propertynames(::HTSRecordFlag)
    propertynames(HTSRecordFlagBits)
end

# TODO: simplify this
function Base.summary(flag::HTSRecordFlag)
    str = [
        flag.paired ? "paired" : "unpaired",
        flag.proper_pair ? "proper_pair" : "improper_pair",
        flag.unmap ? "unmapped" : "mapped",
        flag.munmap ? "mate_unmapped" : "mate_mapped",
        flag.reverse ? "rev_strand" : "pos_strand",
        flag.mreverse ? "mate_rev_strand" : "mate_pos_strand",
        flag.read1 ? "firstread" : "not_firstread",
        flag.read2 ? "lastread" : "not_lastread",
        flag.secondary ? "secondary" : "primary",
        flag.qcfail ? "qcfail" : "not_qcfail",
        flag.dup ? "dup" : "not_dup",
        flag.supplementary ? "supplementary" : "not_supplementary",
    ]
    join(str, ",")
end
