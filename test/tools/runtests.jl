module TestTools

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Seq,
      Bio.Blast,
      TestFunctions

get_bio_fmt_specimens()
path = Pkg.dir("Bio", "test", "BioFmtSpecimens")

@testset "Blast+ blastn" begin
    na1 = dna"CGGACCAGACGGACACAGGGAGAAGCTAGTTTCTTTCATGTGATTGANATNATGACTCTACTCCTAAAAGGGAAAAANCAATATCCTTGTTTACAGAAGAGAAACAAACAAGCCCCACTCAGCTCAGTCACAGGAGAGAN"
    na2 = dna"CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATGGCTGGGAATCTTACTCTTTCATCTGATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATCCTGGTTTTACAGATGTGAAACTTTCAA"
    fna = joinpath(path, "FASTA", "f002.fasta")
    nucldb = joinpath(path, "BLASTDB", "f002")
    nuclresults = joinpath(path, "BLASTDB", "f002.xml")

    @test blastn(na1, na2)
    @test blastn(na1, [na1, na2])
    @test blastn([na1, na2], [na1, na2])
    @test blastn(na1, nucldb, db=true)
    @test blastn(na1, fna)
    @test blastn(fna, nucldb, db=true)

end

@testset "Blast+ blastp" begin
    aa1 = aa"MWATLPLLCAGAWLLGVPVCGAAELSVNSLEKFHFKSWMSKHRKTYSTEEYHHRLQTFASNWRKINAHNNGNHTFKMALNQFSDMSFAEIKHKYLWSEPQNCSATKSNYLRGTGPYPPSVDWRKKGNFVSPVKNQGACGS"
    aa2 = aa"MWTALPLLCAGAWLLSAGATAELTVNAIEKFHFTSWMKQHQKTYSSREYSHRLQVFANNWRKIQAHNQRNHTFKMGLNQFSDMSFAEIKHKYLWSEPQNCSATKSNYLRGTGPYPSSMDWRKKGNVVSPVKNQGACGSCW"
    faa = joinpath(path, "FASTA", "cysprot.fasta")
    protdb = joinpath(path, "BLASTDB", "cysprot")
    protresults = joinpath(path, "BLASTDB", "cysprot.xml")

    @test blastp(aa1, aa2)
    @test blastp(aa1, [aa1, aa2])
    @test blastp([aa1, aa2], [aa1, aa2])
    @test blastp(aa1, protdb, db=true)
    @test blastp(aa1, faa)
    @test blastp(faa, protdb, db=true)
end

end # TestTools
