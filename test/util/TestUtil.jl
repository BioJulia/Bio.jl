module TestUtil

using FactCheck,
    Bio.Util,
    TestFunctions

facts("Sliding-Windows") do
    context("Creation") do
        # Test EachWindow iterators are created correctly when given sensible values.
        # Randomly create some sequences of a random length.
        # Then randomly assign each sequence a winsow size and step size.
        context("Sensible values") do
            for i in 1:100
                n = rand(1:1000)
                testarray = collect(1:n)
                teststring = randstring(n)
                testseq = DNASequence(random_dna(n))
                winsize = rand(1:n)
                stepsize = rand(1:n)
                arrayitr = eachwindow(testarray, winsize, stepsize)
                stritr = eachwindow(teststring, winsize, stepsize)
                seqitr = eachwindow(testseq, winsize, stepsize)
                @fact arrayitr --> EachWindowIterator(testarray, winsize, stepsize)
                @fact stritr --> EachWindowIterator(teststring, winsize, stepsize)
                @fact seqitr --> EachWindowIterator(testseq, winsize, stepsize)
                @fact arrayitr.width --> winsize
                @fact stritr.width --> winsize
                @fact seqitr.width --> winsize
                @fact arrayitr.step --> stepsize
                @fact stritr.step --> stepsize
                @fact seqitr.step --> stepsize
                @fact arrayitr.data --> testarray
                @fact stritr.data --> teststring
                @fact seqitr.data --> testseq
            end
        end
        context("Bad values") do
            for i in 1:100
                n = rand(1:1000)
                testarray = collect(1:n)
                teststring = randstring(n)
                testseq = DNASequence(random_dna(n))
                winsize = rand(n:n + rand(1:1000))
                stepsize = rand(n:n + rand(1:1000))
                @fact_throws eachwindow(testarray, winsize, stepsize)
                @fact_throws eachwindow(teststring, winsize, stepsize)
                @fact_throws eachwindow(testseq, winsize, stepsize)
            end
        end
    end
    context("Iteration") do
        context("Number and size of Windows, and missed elements") do
            for i in 1:100
                n = rand(1:1000)
                testarray = collect(1:n)
                teststring = randstring(n)
                testseq = DNASequence(random_dna(n))
                winsize = rand(1:n)
                stepsize = rand(1:n)
                arrayitr = eachwindow(testarray, winsize, stepsize)
                stringitr = eachwindow(teststring, winsize, stepsize)
                seqitr = eachwindow(testseq, winsize, stepsize)
                arrayres = collect(arrayitr)
                #stringres = collect(stringitr)
                seqres = collect(seqitr)
                expectedMissed = n - StepRange(winsize, stepsize, n)
                @fact length(arrayres) --> size(arrayitr)
                @fact length(stringres) --> size(stringitr)
                @fact length(seqres) --> size(seqitr)
                @fact missed(arrayres) --> expectedMissed
                @fact missed(stringres) --> expectedMissed
                @fact missed(seqres) --> expectedMissed
                for win in arrayres
                    @fact length(win) --> winsize
                end
                #for win in stringres
                #    @fact length(win) --> winsize
                #end
                for win in seqres
                    @fact length(win) --> winsize
                end
            end
        end
    end
end


end
