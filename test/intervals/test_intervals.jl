module TestIntervals

using FactCheck
using Distributions
using Bio
using Bio.Intervals
using YAML
import ..get_bio_fmt_specimens


# Test that an array of intervals is well ordered
function Intervals.isordered{I <: Interval}(intervals::Vector{I})
    for i = 2:length(intervals)
        if !Intervals.isordered(intervals[i-1], intervals[i])
            return false
        end
    end
    return true
end


# Generate an array of n random Interval{Int} object. With sequence names
# samples from seqnames, and intervals drawn to lie in [1, maxpos].
function random_intervals(seqnames, maxpos::Int, n::Int)
    seq_dist = Categorical(length(seqnames))
    strand_dist = Categorical(2)
    length_dist = Normal(1000, 1000)

    intervals = Array(Interval{Int}, n)
    for i in 1:n
        intlen = maxpos
        while intlen >= maxpos || intlen <= 0
            intlen = ceil(Int, rand(length_dist))
        end
        first = rand(1:maxpos-intlen)
        last = first + intlen - 1
        strand = rand(strand_dist) == 1 ? STRAND_POS : STRAND_NEG
        intervals[i] = Interval{Int}(seqnames[rand(seq_dist)],
                                     first, last, strand, i)
    end
    return intervals
end


# A simple interval intersection implementation to test against.
function simple_intersection(intervals_a, intervals_b)
    sort!(intervals_a)
    sort!(intervals_b)

    intersections = Any[]

    i = 1
    j = 1
    while i <= length(intervals_a) && j <= length(intervals_b)
        ai = intervals_a[i]
        bj = intervals_b[j]

        if Intervals.alphanum_isless(ai.seqname, bj.seqname) ||
           (ai.seqname == bj.seqname && ai.last < bj.first)
            i += 1
        elseif Intervals.alphanum_isless(bj.seqname, ai.seqname) ||
               (ai.seqname == bj.seqname && bj.last < ai.first)
            j += 1
        else
            k = j
            while k <= length(intervals_b) && intervals_b[k].first <= ai.last
                if isoverlapping(ai, intervals_b[k])
                    push!(intersections, (ai, intervals_b[k]))
                end
                k += 1
            end
            i += 1
        end
    end

    return intersections
end


function simple_coverage(intervals)
    seqlens = Dict{String, Int}()
    for interval in intervals
        if get(seqlens, interval.seqname, -1) < interval.last
            seqlens[interval.seqname] = interval.last
        end
    end

    covarrays = Dict{String, Vector{Int}}()
    for (seqname, seqlen) in seqlens
        covarrays[seqname] = zeros(Int, seqlen)
    end

    for interval in intervals
        arr = covarrays[interval.seqname]
        for i in interval.first:interval.last
            arr[i] += 1
        end
    end

    covintervals = Interval{Uint32}[]
    for (seqname, arr) in covarrays
        i = j = 1
        while i <= length(arr)
            if arr[i] > 0
                j = i + 1
                while j <= length(arr) && arr[j] == arr[i]
                    j += 1
                end
                push!(covintervals,
                      Interval{Uint32}(seqname, i, j - 1, STRAND_BOTH, arr[i]))
                i = j
            else
                i += 1
            end
        end
    end

    return covintervals
end


facts("IntervalCollection") do

    context("Insertion/Iteration") do
        n = 100000
        intervals = random_intervals(["one", "two", "three"], 1000000, n)
        ic = IntervalCollection{Int}()

        @fact isempty(ic) --> true
        @fact (collect(Interval{Int}, ic) == Interval{Int}[]) --> true

        for interval in intervals
            push!(ic, interval)
        end
        @fact Intervals.isordered(collect(Interval{Int}, ic)) --> true
    end


    context("Intersection") do
        n = 1000
        srand(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "three", "four"], 1000000, n)

        # empty versus empty
        ic_a = IntervalCollection{Int}()
        ic_b = IntervalCollection{Int}()
        @fact (collect(intersect(ic_a, ic_b)) == Any[]) --> true

        # empty versus non-empty
        for interval in intervals_a
            push!(ic_a, interval)
        end

        @fact (collect(intersect(ic_a, ic_b)) == Any[]) --> true
        @fact (collect(intersect(ic_b, ic_a)) == Any[]) --> true

        # non-empty versus non-empty
        for interval in intervals_b
            push!(ic_b, interval)
        end

        @fact (sort(collect(intersect(ic_a, ic_b))) ==
               sort(simple_intersection(intervals_a, intervals_b))) --> true
    end


    context("Show") do
        nullout = open("/dev/null", "w")

        ic = IntervalCollection{Int}()
        show(nullout, ic)

        push!(ic, Interval{Int}("one", 1, 1000, STRAND_POS, 0))
        show(nullout, ic)

        intervals = random_intervals(["one", "two", "three"], 1000000, 100)
        for interval in intervals
            push!(ic, interval)
        end
        show(nullout, ic)

        show(nullout, STRAND_NA)
        show(nullout, STRAND_POS)
        show(nullout, STRAND_NEG)
        show(nullout, STRAND_BOTH)
    end
end


facts("Alphanumeric Sorting") do
    @fact sort(["b", "c" ,"a"], lt=Intervals.alphanum_isless) --> ["a", "b", "c"]
    @fact sort(["a10", "a2" ,"a1"], lt=Intervals.alphanum_isless) --> ["a1", "a2", "a10"]
    @fact sort(["a10a", "a2c" ,"a3b"], lt=Intervals.alphanum_isless) --> ["a2c", "a3b", "a10a"]
    @fact sort(["a3c", "a3b" ,"a3a"], lt=Intervals.alphanum_isless) --> ["a3a", "a3b", "a3c"]
    @fact sort(["a1ac", "a1aa" ,"a1ab"], lt=Intervals.alphanum_isless) --> ["a1aa", "a1ab", "a1ac"]

    @fact Intervals.alphanum_isless("aa", "aa1") --> true
    @fact Intervals.alphanum_isless("aa1", "aa") --> false

    @fact sort(["ac3", "aa", "ab", "aa1", "ac", "ab2"], lt=Intervals.alphanum_isless) -->
        ["aa", "aa1", "ab", "ab2", "ac", "ac3"]
end


facts("IntervalStream") do
    context("StreamBuffer") do
        ref = Int[]
        sb = Intervals.StreamBuffer{Int}()
        @fact isempty(sb) --> true
        @fact (length(sb) == 0) --> true
        @fact_throws shift!(sb)

        ref_shifts = Int[]
        sb_shifts = Int[]

        for i in 1:10000
            if !isempty(ref) && rand() < 0.3
                push!(ref_shifts, shift!(ref))
                push!(sb_shifts, shift!(sb))
            else
                x = rand(Int)
                push!(ref, x)
                push!(sb, x)
            end
        end

        @fact (length(sb) == length(ref)) --> true
        @fact ([sb[i] for i in 1:length(sb)] == ref) --> true
        @fact (ref_shifts == sb_shifts) --> true
        @fact_throws sb[0]
        @fact_throws sb[length(sb) + 1]
    end

    context("Intersection") do
        n = 1000
        srand(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "three", "four"], 1000000, n)

        ic_a = IntervalCollection{Int}()
        ic_b = IntervalCollection{Int}()

        for interval in intervals_a
            push!(ic_a, interval)
        end

        for interval in intervals_b
            push!(ic_b, interval)
        end

        # non-empty versus non-empty, stream intersection
        it = Intervals.IntervalStreamIntersectIterator{Int, Int,
                IntervalCollection{Int}, IntervalCollection{Int}}(
                ic_a, ic_b, Intervals.alphanum_isless)

        @fact (sort(collect(it)) ==
               sort(simple_intersection(intervals_a, intervals_b))) --> true

        # Interesction edge cases: skipping over whole sequences
        typealias SimpleIntersectIterator
            Intervals.IntervalStreamIntersectIterator{Nothing, Nothing,
                Vector{Interval{Nothing}}, Vector{Interval{Nothing}}}

        it = SimpleIntersectIterator(
            [Interval("a", 1, 100, STRAND_POS, nothing), Interval("c", 1, 100, STRAND_POS, nothing)],
            [Interval("a", 1, 100, STRAND_POS, nothing), Interval("b", 1, 100, STRAND_POS, nothing)],
            isless)
        @fact length(collect(it)) --> 1

        it = SimpleIntersectIterator(
            [Interval("c", 1, 100, STRAND_POS, nothing), Interval("d", 1, 100, STRAND_POS, nothing)],
            [Interval("b", 1, 100, STRAND_POS, nothing), Interval("d", 1, 100, STRAND_POS, nothing)],
            isless)
        @fact length(collect(it)) --> 1

        # unsorted streams are not allowed
        @fact_throws begin
            it = SimpleIntersectIterator(
                [Interval("b", 1, 1000, STRAND_POS, nothing),
                 Interval("a", 1, 1000, STRAND_POS, nothing)],
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("b", 1, 1000, STRAND_POS, nothing)], isless)
            collect(it)
        end

        @fact_throws begin
            it = SimpleIntersectIterator(
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("a", 500, 1000, STRAND_POS, nothing),
                 Interval("a", 400, 2000, STRAND_POS, nothing)],
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("b", 1, 1000, STRAND_POS, nothing)], isless)
            collect(it)
        end
    end


    context("IntervalStream Intersection") do
        n = 1000
        srand(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "two", "three"], 1000000, n)

        ic_a = IntervalCollection{Int}()
        ic_b = IntervalCollection{Int}()

        for interval in intervals_a
            push!(ic_a, interval)
        end

        for interval in intervals_b
            push!(ic_b, interval)
        end

        ItType = Intervals.IntervalStreamIntersectIterator{Int, Int,
            IntervalCollection{Int}, IntervalCollection{Int}}

        @fact (sort(collect(ItType(ic_a, ic_b, isless))) ==
               sort(simple_intersection(intervals_a, intervals_b))) --> true
    end

    context("IntervalStream Coverage") do
        n = 10000
        srand(1234)
        intervals = random_intervals(["one", "two", "three"], 1000000, n)

        ic = IntervalCollection{Int}()
        for interval in intervals
            push!(ic, interval)
        end

        @fact (sort(simple_coverage(intervals)) == sort(collect(coverage(ic)))) --> true
    end
end


facts("Interval Parsing") do
    context("BED Parsing") do
        get_bio_fmt_specimens()

        function check_bed_parse(filename)
            # Reading from a stream
            for seqrec in open(open(filename), BED)
            end

            # Reading from a memory mapped file
            for seqrec in open(filename, BED, memory_map=true)
            end

            # Reading from a regular file
            for seqrec in open(filename, BED, memory_map=false)
            end

            return true
        end

        path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "BED")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            valid = get(specimen, "valid", true)
            if valid
                @fact check_bed_parse(joinpath(path, specimen["filename"])) --> true
            else
                @fact_throws check_bed_parse(joinpath(path, specimen["filename"]))
            end
        end
    end

    context("BED Intersection") do
        # Testing strategy: there are two entirely separate intersection
        # algorithms for IntervalCollection and IntervalStream. Here we test
        # them both by checking that they agree by generating and intersecting
        # random BED files.

        function check_intersection(filename_a, filename_b)
            ic_a = IntervalCollection{BEDMetadata}()
            for interval in open(filename_a, BED)
                push!(ic_a, interval)
            end

            ic_b = IntervalCollection{BEDMetadata}()
            for interval in open(filename_b, BED)
                push!(ic_b, interval)
            end

            xs = sort(collect(intersect(open(filename_a, BED), open(filename_b, BED))))
            ys = sort(collect(intersect(ic_a, ic_b)))

            return xs == ys
        end

        n = 10000
        srand(1234)
        intervals_a = random_intervals(["one", "two", "three", "four", "five"], 1000000, n)
        filename_a = Pkg.dir("Bio", "test", "intervals", "test_a.bed")
        out = open(filename_a, "w")
        for interval in sort(intervals_a)
            println(out, interval.seqname, "\t", interval.first - 1, "\t",
                    interval.last, "\t", interval.metadata, "\t", 1000, "\t", interval.strand)
        end
        close(out)

        intervals_b = random_intervals(["one", "two", "three", "four", "five"], 1000000, n)
        filename_b = Pkg.dir("Bio", "test", "intervals", "test_b.bed")
        ic_b = IntervalCollection{Int}()
        out = open(filename_b, "w")
        for interval in sort(intervals_b)
            println(out, interval.seqname, "\t", interval.first - 1, "\t",
                    interval.last, "\t", interval.metadata, "\t", 1000, "\t", interval.strand)
            push!(ic_b, interval)
        end
        close(out)

        @fact check_intersection(filename_a, filename_b) --> true
    end
end


#facts("BigBed") do
    #context("BED → BigBed → BED round-trip") do
        #path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "BED")
        #for specimen in YAML.load_file(joinpath(path, "index.yml"))
            #if !get(specimen, "valid", true)
                #continue
            #end

            ## BED → BigBed
            #intervals = IntervalCollection(
                #read(open(joinpath(path, specimen["filename"])), BED))
            #out = IOBuffer()
            #write(out, BigBed, intervals)
            #bigbed_data = takebuf_array(out)

            ## BigBed → BED
            #bb = read(bigbed_data, BigBed)
            #intervals2 = IntervalCollection(bb)

            #@fact intervals == intervals2 --> true
        #end
    #end

    #context("BigBed Intersection") do
        #n = 10000
        #srand(1234)
        #chroms = ["one", "two", "three", "four", "five"]
        #intervals = IntervalCollection(
            #[Interval{BEDMetadata}(i.seqname, i.first, i.last, i.strand, BEDMetadata("", 1000))
             #for i in random_intervals(chroms, 1000000, n)], true)

        ## convert to bigbed in memory
        #out = IOBuffer()
        #write(out, BigBed, intervals)
        #bb = read(takebuf_array(out), BigBed)

        ## intersection queries
        #num_queries = 1000
        #queries = random_intervals(chroms, 1000000, num_queries)

        #@fact all(Bool[IntervalCollection(collect(BEDInterval, intersect(intervals, query))) ==
                       #IntervalCollection(collect(BEDInterval, intersect(bb, query)))
                       #for query in queries]) --> true
    #end

    ## TODO: test summary information against output from kent's bigBedSummary
#end

end # module TestIntervals
