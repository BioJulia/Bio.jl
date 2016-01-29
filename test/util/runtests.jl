module TestUtil

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Util,
    Bio.Seq,
    TestFunctions

function generate_testcase_data()
    n = rand(1:1000)
    arr = collect(1:n)
    str = randstring(n)
    nuc = DNASequence(random_dna(n))
    rib = RNASequence(random_rna(n))
    aas = AminoAcidSequence(random_aa(n))
    return n, arr, str, nuc, rib, aas
end

function generate_test_window(n::Int)
    return rand(1:n), rand(1:n)
end

function generate_test_itrs(width, step, arr, str, nuc, rib, aas)
    arritr = eachwindow(arr, width, step)
    stritr = eachwindow(str, width, step)
    nucitr = eachwindow(nuc, width, step)
    ribitr = eachwindow(rib, width, step)
    aasitr = eachwindow(aas, width, step)
    return arritr, stritr, nucitr, ribitr, aasitr
end

@testset "Sliding-Windows" begin
    @testset "Iterator Creation" begin
        @testset "Sensible values" begin
            for _ in 1:10
                n, arr, str, nuc, rib, aas = generate_testcase_data()
                winwidth, winstep = generate_test_window(n)
                arritr, stritr, nucitr, ribitr, aasitr =
                    generate_test_itrs(winwidth, winstep, arr, str, nuc, rib, aas)
                @test arritr == EachWindowIterator{Array{Int, 1}}(arr, winwidth, winstep)
                @test stritr == EachWindowIterator{ASCIIString}(str, winwidth, winstep)
                @test nucitr == EachWindowIterator{DNASequence}(nuc, winwidth, winstep)
                @test ribitr == EachWindowIterator{RNASequence}(rib, winwidth, winstep)
                @test aasitr == EachWindowIterator{AminoAcidSequence}(aas, winwidth, winstep)
                @test arritr.width == winwidth
                @test stritr.width == winwidth
                @test nucitr.width == winwidth
                @test ribitr.width == winwidth
                @test aasitr.width == winwidth
                @test arritr.step == winstep
                @test stritr.step == winstep
                @test nucitr.step == winstep
                @test ribitr.step == winstep
                @test aasitr.step == winstep
                @test arritr.data == arr
                @test stritr.data == str
                @test nucitr.data == nuc
                @test ribitr.data == rib
                @test aasitr.data == aas
            end
        end

        @testset "Bad values" begin
            for _ in 1:10
                n, arr, str, nuc = generate_testcase_data()
                winwidth = rand(n + 1:n + rand(1:1000))
                winstep = rand(n + 1:n + rand(1:1000))
                @test_throws Exception eachwindow(arr, winwidth, winstep)
                @test_throws Exception eachwindow(str, winwidth, winstep)
                @test_throws Exception eachwindow(nuc, winwidth, winstep)
                @test_throws Exception eachwindow(rib, winwidth, winstep)
                @test_throws Exception eachwindow(aas, winwidth, winstep)
            end
        end
    end


    @testset "Iteration" begin
        @testset "N & size of Windows, N of missed elements" begin
            for _ in 1:10
                n, arr, str, nuc, rib, aas = generate_testcase_data()
                winwidth, winstep = generate_test_window(n)
                arritr, stritr, nucitr, ribitr, aasitr =
                    generate_test_itrs(winwidth, winstep, arr, str, nuc, rib, aas)
                arrres = collect(arritr)
                strres = collect(stritr)
                nucres = collect(nucitr)
                ribres = collect(ribitr)
                aasres = collect(aasitr)
                expectedMissed = n - StepRange(winwidth, winstep, n).stop
                @test length(arrres) == size(arritr)
                @test length(strres) == size(stritr)
                @test length(nucres) == size(nucitr)
                @test length(ribres) == size(ribitr)
                @test length(aasres) == size(aasitr)
                @test missed(arritr) == expectedMissed
                @test missed(stritr) == expectedMissed
                @test missed(nucitr) == expectedMissed
                @test missed(ribitr) == expectedMissed
                @test missed(aasitr) == expectedMissed
                for win in arrres
                    @test length(win) == winwidth
                end
                for win in strres
                    @test length(win) == winwidth
                end
                for win in nucres
                    @test length(win) == winwidth
                end
                for win in ribres
                    @test length(win) == winwidth
                end
                for win in aasres
                    @test length(win) == winwidth
                end
            end
        end

        @testset "Content of windows generated" begin
            @testset "DNA" begin
                nuc = DNASequence(random_dna(10))
                for width = 1:5, step = 1:5
                    nucitr = eachwindow(nuc, width, step)
                    i, j = 1, width
                    for win in nucitr
                        @test win == sub(nuc, i:j)
                        @test win == nuc[i:j]
                        i += step
                        j += step
                    end
                end
            end

            @testset "RNA" begin
                rib = RNASequence(random_rna(10))
                for width = 1:5, step = 1:5
                    ribitr = eachwindow(rib, width, step)
                    i, j = 1, width
                    for win in ribitr
                        @test win == sub(rib, i:j)
                        @test win == rib[i:j]
                        i += step
                        j += step
                    end
                end
            end

            @testset "Amino Acids" begin
                aas = AminoAcidSequence(random_aa(10))
                for width = 1:5, step = 1:5
                    aasitr = eachwindow(aas, width, step)
                    i, j = 1, width
                    for win in aasitr
                        @test win == sub(aas, i:j)
                        @test win == aas[i:j]
                        i += step
                        j += step
                    end
                end
            end

            @testset "Arrays" begin
                arr = collect(1:10)
                for width = 1:5, step = 1:5
                    arritr = eachwindow(arr, width, step)
                    i, j = 1, width
                    for win in arritr
                        @test win == sub(arr, i:j)
                        i += step
                        j += step
                    end
                end
            end

            @testset "Strings" begin
                str = randstring(10)
                for width = 1:5, step = 1:5
                    stritr = eachwindow(str, width, step)
                    i, j = 1, width
                    for win in stritr
                        @test win == str[i:j]
                        i += step
                        j += step
                    end
                end
            end

        end
    end
end

end
