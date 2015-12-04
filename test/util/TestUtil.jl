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
                seq = collect(1:n)
                winsize = rand(1:length(seq))
                stepsize = rand(1:length(seq))
                itr = eachwindow(seq, winsize, stepsize)
                @fact itr --> EachWindowIterator(seq, winsize, stepsize)
                @fact itr.width --> winsize
                @fact itr.step --> stepsize
                @fact itr.data --> seq
            end
        end
        context("Bad values") do
            for i in 1:100
                n = rand(1:1000)
                seq = collect(1:n)
                winsize = rand(length(seq):length(seq) + rand(1:1000))
                stepsize = rand(length(seq):length(seq) + rand(1:1000))
                @fact_throws eachwindow(seq, winsize, stepsize)
            end
        end
    end
    context("Iteration") do
        context("Number and size of Windows") do
            for i in 1:100
                n = rand(1:1000)
                seq = collect(1:n)
                winsize = rand(1:length(seq))
                stepsize = rand(1:length(seq))
                itr = eachwindow(seq, winsize, stepsize)
                res = collect(itr)
                @fact length(res) --> size(itr)
                for win in res
                    @fact length(res) --> winsize
                end
            end
        end
    end
end


end
