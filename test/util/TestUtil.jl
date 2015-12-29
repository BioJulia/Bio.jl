module TestUtil

using FactCheck,
    Bio.Util,
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

facts("Sliding-Windows") do
    context("Iterator Creation") do
        context("Sensible values") do
            for _ in 1:10
                n, arr, str, nuc, rib, aas = generate_testcase_data()
                winwidth, winstep = generate_test_window(n)
                arritr, stritr, nucitr, ribitr, aasitr =
                    generate_test_itrs(winwidth, winstep, arr, str, nuc, rib, aas)
                @fact arritr --> EachWindowIterator{Array{Int, 1}}(arr, winwidth, winstep)
                @fact stritr --> EachWindowIterator{ASCIIString}(str, winwidth, winstep)
                @fact nucitr --> EachWindowIterator{DNASequence}(nuc, winwidth, winstep)
                @fact ribitr --> EachWindowIterator{RNASequence}(rib, winwidth, winstep)
                @fact aasitr --> EachWindowIterator{AminoAcidSequence}(aas, winwidth, winstep)
                @fact arritr.width --> winwidth
                @fact stritr.width --> winwidth
                @fact nucitr.width --> winwidth
                @fact ribitr.width --> winwidth
                @fact aasitr.width --> winwidth
                @fact arritr.step --> winstep
                @fact stritr.step --> winstep
                @fact nucitr.step --> winstep
                @fact ribitr.step --> winstep
                @fact aasitr.step --> winstep
                @fact arritr.data --> arr
                @fact stritr.data --> str
                @fact nucitr.data --> nuc
                @fact ribitr.data --> rib
                @fact aasitr.data --> aas
            end
        end

        context("Bad values") do
            for _ in 1:10
                n, arr, str, nuc = generate_testcase_data()
                winwidth = rand(n + 1:n + rand(1:1000))
                winstep = rand(n + 1:n + rand(1:1000))
                @fact_throws eachwindow(arr, winwidth, winstep)
                @fact_throws eachwindow(str, winwidth, winstep)
                @fact_throws eachwindow(nuc, winwidth, winstep)
                @fact_throws eachwindow(rib, winwidth, winstep)
                @fact_throws eachwindow(aas, winwidth, winstep)
            end
        end
    end


    context("Iteration") do
        context("N & size of Windows, N of missed elements") do
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
                @fact length(arrres) --> size(arritr)
                @fact length(strres) --> size(stritr)
                @fact length(nucres) --> size(nucitr)
                @fact length(ribres) --> size(ribitr)
                @fact length(aasres) --> size(aasitr)
                @fact missed(arritr) --> expectedMissed
                @fact missed(stritr) --> expectedMissed
                @fact missed(nucitr) --> expectedMissed
                @fact missed(ribitr) --> expectedMissed
                @fact missed(aasitr) --> expectedMissed
                for win in arrres
                    @fact length(win) --> winwidth
                end
                for win in strres
                    @fact length(win) --> winwidth
                end
                for win in nucres
                    @fact length(win) --> winwidth
                end
                for win in ribres
                    @fact length(win) --> winwidth
                end
                for win in aasres
                    @fact length(win) --> winwidth
                end
            end
        end

        context("Content of windows generated") do
            context("DNA") do
                nuc = DNASequence(random_dna(10))
                for width = 1:5, step = 1:5
                    nucitr = eachwindow(nuc, width, step)
                    i, j = 1, width
                    for win in nucitr
                        @fact win --> sub(nuc, i:j)
                        @fact win --> nuc[i:j]
                        i += step
                        j += step
                    end
                end
            end

            context("RNA") do
                rib = RNASequence(random_rna(10))
                for width = 1:5, step = 1:5
                    ribitr = eachwindow(rib, width, step)
                    i, j = 1, width
                    for win in ribitr
                        @fact win --> sub(rib, i:j)
                        @fact win --> rib[i:j]
                        i += step
                        j += step
                    end
                end
            end

            context("Amino Acids") do
                aas = AminoAcidSequence(random_aa(10))
                for width = 1:5, step = 1:5
                    aasitr = eachwindow(aas, width, step)
                    i, j = 1, width
                    for win in aasitr
                        @fact win --> sub(aas, i:j)
                        @fact win --> aas[i:j]
                        i += step
                        j += step
                    end
                end
            end

            context("Arrays") do
                arr = collect(1:10)
                for width = 1:5, step = 1:5
                    arritr = eachwindow(arr, width, step)
                    i, j = 1, width
                    for win in arritr
                        @fact win --> sub(arr, i:j)
                        i += step
                        j += step
                    end
                end
            end

            context("Strings") do
                str = randstring(10)
                for width = 1:5, step = 1:5
                    stritr = eachwindow(str, width, step)
                    i, j = 1, width
                    for win in stritr
                        @fact win --> str[i:j]
                        i += step
                        j += step
                    end
                end
            end

        end
    end
end

end
