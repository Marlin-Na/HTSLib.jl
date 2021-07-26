using HTSLib
using HTSLib: htslib  # low-level API module
using Test
using FormatSpecimens
import HttpIO

# Low-level API
# -------------

@testset "Low-level API" begin
    # check struct size (these sizes are derived from C's sizeof operator)
    if Sys.isunix() && Sys.ARCH == :x86_64  # LP64
        @test sizeof(htslib.kstring_t) == 24
        @test sizeof(htslib.htsFormat) == 32
        @test sizeof(htslib.htsFile) == 128
        @test sizeof(htslib.bam1_core_t) == 48
        @test sizeof(htslib.bam1_t) == 80
    end
end

# High-level API
# --------------

@testset "HTSLib version" begin
    @test HTSLib.HTSLIB_VERSION == v"1.10.2"
end

@testset "SAM header" begin
    bamfile = joinpath(path_of_format("BAM"), "ce#1.bam")
    reader = open(HTSReadWriter, bamfile)
    header_str = string(header(reader))

    @test isa(header_str, String)
    @test string(SamHeader(header_str)) == header_str
end

@testset "BAM reading" begin
    @testset "Read local bam" begin
        bamfile = joinpath(path_of_format("BAM"), "ce#1.bam")
        reader = open(HTSReadWriter, bamfile)

        @test eltype(reader) == HTSRecord
        @test isa(header(reader), SamHeader)

        record = HTSRecord()
        read!(reader, record)

        @test refid(record) == 1
        @test leftposition(record) == 2
        @test rightposition(record) == 102
        @test sequence(record, 3) == 'T'
        @test seqlength(record) == 100
        @test sequence(record) == collect("CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA")
        @test quality(record) == [Int(x) - 33 for x in "#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"]
        @test eltype(quality(record)) == UInt8
    end

    @testset "iterator" begin
        bamfile1 = joinpath(path_of_format("BAM"), "ce#1.bam")
        bamfile2 = joinpath(path_of_format("BAM"), "ce#2.bam")

        @test length(collect(HTSReadWriter(bamfile1))) == 1
        @test length(collect(HTSReadWriter(bamfile2))) == 2
    end
end

@testset "BAM reading from url" begin
    ## TODO: test invalid url
    bamfile1 = "https://raw.githubusercontent.com/BioJulia/FormatSpecimens.jl/master/BAM/ce#1.bam"
    bamfile2 = "https://raw.githubusercontent.com/BioJulia/FormatSpecimens.jl/master/BAM/ce#2.bam"

    @test length(collect(HTSReadWriter(bamfile1))) == 1
    @test length(collect(HTSReadWriter(bamfile2))) == 2
end

@testset "BAM/CRAM reading from gsurl" begin
    # skip test if gcloud is not installed
    if isnothing(Sys.which("gcloud"))
        return
    end

    gcsclient = HttpIO.GCSClient()

    cram = "gs://gatk-test-data/wgs_cram/NA12878_20k_hg38/NA12878.cram"
    crai = "gs://gatk-test-data/wgs_cram/NA12878_20k_hg38/NA12878.crai"

    bam = "gs://gatk-test-data/wgs_bam/NA12878_20k_hg38/NA12878.bam"
    bai = "gs://gatk-test-data/wgs_bam/NA12878_20k_hg38/NA12878.bai"

    @testset "Simple BAM read" begin
        bamio = HTSReadWriter(bam; index=bai, gcsclient=gcsclient)
        record = HTSRecord()
        read!(bamio, record)
        @test leftposition(record) == 10001
        @test rightposition(record) == 10110
    end

    @testset "Simple CRAM read" begin
        cramio = HTSReadWriter(cram; gcsclient=gcsclient) # TODO: hts_idx_load2 does not recognize cram index file
        record = HTSRecord()
        read!(cramio, record)
        @test leftposition(record) == 10001
        @test rightposition(record) == 10110
    end

    @testset "Bam random access" begin
        bamio = HTSReadWriter(bam; index=bai, gcsclient=gcsclient)

        @testset "single region" begin
            bamitr1 = HTSRegionsIterator(bamio, "chr1:1000000-2000000")
            bamitr2 = HTSRegionsIterator(bamio, "chr1", 1000000, 2000000)

            @test length(collect(bamitr1)) == 25
            @test length(collect(bamitr2)) == 25

            @test !isnothing(tryread!(bamio, HTSRecord()))
        end

        @testset "unmapped reads" begin
            bamitr = HTSRegionsIterator(bamio, "*")
            @test length(collect(bamitr)) == 442

            ## we have iterated unmapped reads at end of the file
            @test isnothing(tryread!(bamio, HTSRecord()))
        end

        @testset "GC" begin
            itr = open(HTSReadWriter, bam; index=bai, gcsclient=gcsclient) do reader
                itr = HTSRegionsIterator(reader, "chr1", 1000000, 2000000)
                @test !isnothing(read(itr))
                itr
            end
            # reader has been closed, itr is invalid now, but it should be safe to finalize it.
            finalize(itr)
        end
    end

    @testset "CRAM random access" begin
        # TODO: add index
        #cramio = HTSReadWriter(cram; gcsclient=gcsclient)

        # @testset "from file start" begin
        #     itr = open(HTSReadWriter, cram; gcsclient=gcsclient) do reader
        #         itr = HTSRegionsIterator(reader, ".")
        #         @test !isnothing(read(itr))
        #         itr
        #     end
        #     finalize(itr)
        # end
    end
end

@testset "Write hts file" begin
    @testset "bam" begin
        bamfile = joinpath(path_of_format("BAM"), "ce#1.bam")
        bamin = open(HTSReadWriter, bamfile)
        (path, io) = mktemp()
        bamout = HTSReadWriter(io, "wbz")

        # write header
        write(bamout, header(bamin))
        # write records
        for r in bamin
            write(bamout, r)
        end

        # flush and close
        close(bamout)

        # reread
        f = HTSReadWriter(path)
        @test HTSLib.format_description(f) == "BAM version 1 compressed sequence data"
        records = collect(f)
        @test length(records) == 1
        @test refid(records[1]) == 1
        @test rightposition(records[1]) == 102
    end
end

@testset "Regions iterator" begin

    @testset "Constructor" begin
        
        bam = joinpath(path_of_format("BAM"), "R_12h_D06.uniq.q40.bam")
        bai = joinpath(path_of_format("BAM"), "R_12h_D06.uniq.q40.bam.bai")
        reader = HTSReadWriter(bam; index=bai)
        
        # Vector{String}
        @test length(collect(HTSRegionsIterator(reader, ["chrM:1-500"]))) == 2
        # String
        @test length(collect(HTSRegionsIterator(reader, "chrM:1-500"))) == 2
        # Tuple
        @test length(collect(HTSRegionsIterator(reader, ("chrM", 1, 500)))) == 2
        @test length(collect(HTSRegionsIterator(reader, (22, 1, 500)))) == 2
        # Vector{Tuple}
        @test length(collect(HTSRegionsIterator(reader, [("chrM", 1, 500)]))) == 2
        @test length(collect(HTSRegionsIterator(reader, [(22, 1, 500)]))) == 2
        # args
        @test length(collect(HTSRegionsIterator(reader, "chrM", 1, 500))) == 2
        @test length(collect(HTSRegionsIterator(reader, 22, 1, 500))) == 2
        
    end
end
