module TestIntervals

using FactCheck
using Distributions
using Bio
using Bio.Intervals


function random_intervals(seqnames, maxpos::Int, n::Int)
    seq_dist = Categorical(length(seqnames))
    strand_dist = Categorical(2)
    length_dist = Normal(1000, 1000)

    intervals = Array(Interval{Int}, n)
    for i in 1:n
        intlen = ceil(Int, rand(length_dist))
        first = rand(1:maxpos-intlen)
        last = first + intlen - 1
        strand = rand(strand_dist) == 1 ? STRAND_POS : STRAND_NEG
        intervals[i] = Interval{Int}(seqnames[rand(seq_dist)],
                                     first, last, strand, i)
    end
    return intervals
end

facts("IntervalSet Insertion") do
    n = 100000
    intervals = random_intervals(["one", "two", "three"], 1000000, n)
    is = IntervalSet{Int}()
    for interval in intervals
        push!(is, interval)
    end
    @fact sort(collect(is)) == sort(intervals) => true
end

end


