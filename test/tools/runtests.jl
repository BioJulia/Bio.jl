module TestTools

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Compat
using Bio.Seq,
      Bio.Tools.BLAST,
      TestFunctions

get_bio_fmt_specimens()
path = Pkg.dir("Bio", "test", "BioFmtSpecimens")

if !is_windows()  # temporarily disable the BLAST tests on Windows (issue: #197)

@testset "BLAST+ blastn" begin
    na1 = dna"""
    CGGACCAGACGGACACAGGGAGAAGCTAGTTTCTTTCATGTGATTGANAT
    NATGACTCTACTCCTAAAAGGGAAAAANCAATATCCTTGTTTACAGAAGA
    GAAACAAACAAGCCCCACTCAGCTCAGTCACAGGAGAGAN
    """
    na2 = dna"""
    CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATG
    GCTGGGAATCTTACTCTTTCATCTGATACCTTGTTCAGATTTCAAAATAG
    TTGTAGCCTTATCCTGGTTTTACAGATGTGAAACTTTCAA
    """
    fna = joinpath(path, "FASTA", "f002.fasta")
    nucldb = joinpath(path, "BLASTDB", "f002")
    nuclresults = joinpath(path, "BLASTDB", "f002.xml")

    @test typeof(blastn(na1, na2)) == Array{BLASTResult, 1}
    @test typeof(blastn(na1, [na1, na2])) == Array{BLASTResult, 1}
    @test typeof(blastn([na1, na2], [na1, na2])) == Array{BLASTResult, 1}
    @test typeof(blastn(na1, nucldb, db=true)) == Array{BLASTResult, 1}
    @test typeof(blastn(na1, fna)) == Array{BLASTResult, 1}
    @test typeof(blastn(fna, nucldb, db=true)) == Array{BLASTResult, 1}
end

@testset "BLAST+ blastp" begin
    aa1 = aa"""
    MWATLPLLCAGAWLLGVPVCGAAELSVNSLEKFHFKSWMSKHRKTYSTEE
    YHHRLQTFASNWRKINAHNNGNHTFKMALNQFSDMSFAEIKHKYLWSEPQ
    NCSATKSNYLRGTGPYPPSVDWRKKGNFVSPVKNQGACGS
    """
    aa2 = aa"""
    MWTALPLLCAGAWLLSAGATAELTVNAIEKFHFTSWMKQHQKTYSSREYS
    HRLQVFANNWRKIQAHNQRNHTFKMGLNQFSDMSFAEIKHKYLWSEPQNC
    SATKSNYLRGTGPYPSSMDWRKKGNVVSPVKNQGACGSCW
    """
    faa = joinpath(path, "FASTA", "cysprot.fasta")
    protdb = joinpath(path, "BLASTDB", "cysprot")
    protresults = joinpath(path, "BLASTDB", "cysprot.xml")

    @test typeof(blastp(aa1, aa2)) == Array{BLASTResult, 1}
    @test typeof(blastp(aa1, [aa1, aa2])) == Array{BLASTResult, 1}
    @test typeof(blastp([aa1, aa2], [aa1, aa2])) == Array{BLASTResult, 1}
    @test typeof(blastp(aa1, protdb, db=true)) == Array{BLASTResult, 1}
    @test typeof(blastp(aa1, faa)) == Array{BLASTResult, 1}
    @test typeof(blastp(faa, protdb, db=true)) == Array{BLASTResult, 1}
end

end  # if !is_windows()

end # TestTools
