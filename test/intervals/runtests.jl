module TestIntervals

using Base.Test

using Bio.Intervals,
    Distributions,
    YAML,
    TestFunctions

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
    seqlens = Dict{AbstractString, Int}()
    for interval in intervals
        if get(seqlens, interval.seqname, -1) < interval.last
            seqlens[interval.seqname] = interval.last
        end
    end

    covarrays = Dict{AbstractString, Vector{Int}}()
    for (seqname, seqlen) in seqlens
        covarrays[seqname] = zeros(Int, seqlen)
    end

    for interval in intervals
        arr = covarrays[interval.seqname]
        for i in interval.first:interval.last
            arr[i] += 1
        end
    end

    covintervals = Interval{UInt32}[]
    for (seqname, arr) in covarrays
        i = j = 1
        while i <= length(arr)
            if arr[i] > 0
                j = i + 1
                while j <= length(arr) && arr[j] == arr[i]
                    j += 1
                end
                push!(covintervals,
                      Interval{UInt32}(seqname, i, j - 1, STRAND_BOTH, arr[i]))
                i = j
            else
                i += 1
            end
        end
    end

    return covintervals
end

@testset "Strand" begin
    @testset "Constructor" begin
        @test Strand('?') === STRAND_NA
        @test Strand('+') === STRAND_POS
        @test Strand('-') === STRAND_NEG
        @test Strand('.') === STRAND_BOTH
        @test_throws Exception Strand('x')
    end

    @testset "Conversion" begin
        @test convert(Strand, '?') === STRAND_NA
        @test convert(Strand, '+') === STRAND_POS
        @test convert(Strand, '-') === STRAND_NEG
        @test convert(Strand, '.') === STRAND_BOTH

        @test convert(Char, STRAND_NA) === '?'
        @test convert(Char, STRAND_POS) === '+'
        @test convert(Char, STRAND_NEG) === '-'
        @test convert(Char, STRAND_BOTH) === '.'
    end

    @testset "Order" begin
        @test STRAND_NA < STRAND_POS < STRAND_NEG < STRAND_BOTH
    end

    @testset "Show" begin
        @testset "show" begin
            buf = IOBuffer()
            for s in [STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH]
                show(buf, s); print(buf, " ")
            end
            @test takebuf_string(buf) == "STRAND_NA STRAND_POS STRAND_NEG STRAND_BOTH "
        end

        @testset "print" begin
            buf = IOBuffer()
            for s in [STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH]
                print(buf, s)
            end
            @test takebuf_string(buf) == "?+-."
        end
    end
end

@testset "Interval" begin
    @testset "Constructor" begin
        i = Interval("chr1", 10, 20)
        @test seqname(i) == "chr1"
        @test leftposition(i) == 10
        @test rightposition(i) == 20
        @test strand(i) == STRAND_BOTH

        i1 = Interval("chr1", 10, 20, '+')
        i2 = Interval("chr1", 10, 20, STRAND_POS)
        @test i1 == i2

        i1 = Interval("chr2", 5692667, 5701385, '+',        "SOX11")
        i2 = Interval("chr2", 5692667, 5701385, STRAND_POS, "SOX11")
        @test i1 == i2
    end
end

@testset "IntervalCollection" begin

    @testset "Insertion/Iteration" begin
        n = 100000
        intervals = random_intervals(["one", "two", "three"], 1000000, n)
        ic = IntervalCollection{Int}()

        @test isempty(ic)
        @test collect(Interval{Int}, ic) == Interval{Int}[]

        for interval in intervals
            push!(ic, interval)
        end
        @test Intervals.isordered(collect(Interval{Int}, ic))
    end


    @testset "Intersection" begin
        n = 1000
        srand(1234)
        intervals_a = random_intervals(["one", "two", "three"], 1000000, n)
        intervals_b = random_intervals(["one", "three", "four"], 1000000, n)

        # empty versus empty
        ic_a = IntervalCollection{Int}()
        ic_b = IntervalCollection{Int}()
        @test collect(intersect(ic_a, ic_b)) == Any[]

        # empty versus non-empty
        for interval in intervals_a
            push!(ic_a, interval)
        end

        @test collect(intersect(ic_a, ic_b)) == Any[]
        @test collect(intersect(ic_b, ic_a)) == Any[]

        # non-empty versus non-empty
        for interval in intervals_b
            push!(ic_b, interval)
        end

        @test sort(collect(intersect(ic_a, ic_b))) == sort(simple_intersection(intervals_a, intervals_b))
    end


    @testset "Show" begin
        ic = IntervalCollection{Int}()
        show(DevNull, ic)

        push!(ic, Interval{Int}("one", 1, 1000, STRAND_POS, 0))
        show(DevNull, ic)

        intervals = random_intervals(["one", "two", "three"], 1000000, 100)
        for interval in intervals
            push!(ic, interval)
        end
        show(DevNull, ic)

        show(DevNull, STRAND_NA)
        show(DevNull, STRAND_POS)
        show(DevNull, STRAND_NEG)
        show(DevNull, STRAND_BOTH)
    end
end


@testset "Alphanumeric Sorting" begin
    @test sort(["b", "c" ,"a"], lt=Intervals.alphanum_isless) == ["a", "b", "c"]
    @test sort(["a10", "a2" ,"a1"], lt=Intervals.alphanum_isless) == ["a1", "a2", "a10"]
    @test sort(["a10a", "a2c" ,"a3b"], lt=Intervals.alphanum_isless) == ["a2c", "a3b", "a10a"]
    @test sort(["a3c", "a3b" ,"a3a"], lt=Intervals.alphanum_isless) == ["a3a", "a3b", "a3c"]
    @test sort(["a1ac", "a1aa" ,"a1ab"], lt=Intervals.alphanum_isless) == ["a1aa", "a1ab", "a1ac"]

    @test Intervals.alphanum_isless("aa", "aa1")
    @test !Intervals.alphanum_isless("aa1", "aa")

    @test sort(["ac3", "aa", "ab", "aa1", "ac", "ab2"], lt=Intervals.alphanum_isless) ==
        ["aa", "aa1", "ab", "ab2", "ac", "ac3"]
end


@testset "IntervalStream" begin
    @testset "StreamBuffer" begin
        ref = Int[]
        sb = Intervals.StreamBuffer{Int}()
        @test isempty(sb)
        @test length(sb) == 0
        @test_throws Exception shift!(sb)

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

        @test length(sb) == length(ref)
        @test [sb[i] for i in 1:length(sb)] == ref
        @test ref_shifts == sb_shifts
        @test_throws BoundsError sb[0]
        @test_throws BoundsError sb[length(sb) + 1]
    end

    @testset "Intersection" begin
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

        @test sort(collect(it)) == sort(simple_intersection(intervals_a, intervals_b))

        # Interesction edge cases: skipping over whole sequences
        typealias SimpleIntersectIterator
            Intervals.IntervalStreamIntersectIterator{Void, Void,
                Vector{Interval{Void}}, Vector{Interval{Void}}}

        it = SimpleIntersectIterator(
            [Interval("a", 1, 100, STRAND_POS, nothing), Interval("c", 1, 100, STRAND_POS, nothing)],
            [Interval("a", 1, 100, STRAND_POS, nothing), Interval("b", 1, 100, STRAND_POS, nothing)],
            isless)
        @test length(collect(it)) == 1

        it = SimpleIntersectIterator(
            [Interval("c", 1, 100, STRAND_POS, nothing), Interval("d", 1, 100, STRAND_POS, nothing)],
            [Interval("b", 1, 100, STRAND_POS, nothing), Interval("d", 1, 100, STRAND_POS, nothing)],
            isless)
        @test length(collect(it)) == 1

        # unsorted streams are not allowed
        @test_throws Exception begin
            it = SimpleIntersectIterator(
                [Interval("b", 1, 1000, STRAND_POS, nothing),
                 Interval("a", 1, 1000, STRAND_POS, nothing)],
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("b", 1, 1000, STRAND_POS, nothing)], isless)
            collect(it)
        end

        @test_throws Exception begin
            it = SimpleIntersectIterator(
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("a", 500, 1000, STRAND_POS, nothing),
                 Interval("a", 400, 2000, STRAND_POS, nothing)],
                [Interval("a", 1, 1000, STRAND_POS, nothing),
                 Interval("b", 1, 1000, STRAND_POS, nothing)], isless)
            collect(it)
        end
    end


    @testset "IntervalStream Intersection" begin
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

        @test sort(collect(ItType(ic_a, ic_b, isless))) == sort(simple_intersection(intervals_a, intervals_b))
    end

    @testset "IntervalStream Coverage" begin
        n = 10000
        srand(1234)
        intervals = random_intervals(["one", "two", "three"], 1000000, n)

        ic = IntervalCollection{Int}()
        for interval in intervals
            push!(ic, interval)
        end

        @test sort(simple_coverage(intervals)) == sort(collect(coverage(ic)))
    end
end


@testset "Interval Parsing" begin
    @testset "BED Parsing" begin
        get_bio_fmt_specimens()

        function check_bed_parse(filename)
            # Reading from a stream
            for interval in BEDReader(open(filename))
            end

            # Reading from a regular file
            for interval in open(BEDReader, filename)
            end

            # in-place parsing
            stream = open(BEDReader, filename)
            entry = eltype(stream)()
            while !eof(stream)
                read!(stream, entry)
            end
            close(stream)

            # Check round trip
            output = IOBuffer()
            writer = BEDWriter(output)
            expected_entries = Any[]
            for interval in open(BEDReader, filename)
                write(writer, interval)
                push!(expected_entries, interval)
            end
            flush(writer)

            seekstart(output)
            read_entries = Any[]
            for interval in BEDReader(output)
                push!(read_entries, interval)
            end

            return expected_entries == read_entries
        end

        path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "BED")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            valid = get(specimen, "valid", true)
            if valid
                @test check_bed_parse(joinpath(path, specimen["filename"]))
            else
                @test_throws Exception check_bed_parse(joinpath(path, specimen["filename"]))
            end
        end
    end

    @testset "BED Intersection" begin
        # Testing strategy: there are two entirely separate intersection
        # algorithms for IntervalCollection and IntervalStream. Here we test
        # them both by checking that they agree by generating and intersecting
        # random BED files.

        function check_intersection(filename_a, filename_b)
            ic_a = IntervalCollection{BEDMetadata}()
            open(BEDReader, filename_a) do intervals
                for interval in intervals
                    push!(ic_a, interval)
                end
            end

            ic_b = IntervalCollection{BEDMetadata}()
            open(BEDReader, filename_b) do intervals
                for interval in intervals
                    push!(ic_b, interval)
                end
            end

            # This is refactored out to close streams
            fa = open(BEDReader, filename_a)
            fb = open(BEDReader, filename_b)
            xs = sort(collect(intersect(fa, fb)))
            close(fa)
            close(fb)

            ys = sort(collect(intersect(ic_a, ic_b)))

            return xs == ys
        end

        function write_intervals(filename, intervals)
            open(filename, "w") do out
                for interval in sort(intervals)
                    println(out, interval.seqname, "\t", interval.first - 1,
                            "\t", interval.last, "\t", interval.metadata, "\t",
                            1000, "\t", interval.strand)
                end
            end

        end

        n = 10000
        srand(1234)
        intervals_a = random_intervals(["one", "two", "three", "four", "five"],
                                       1000000, n)
        intervals_b = random_intervals(["one", "two", "three", "four", "five"],
                                       1000000, n)

        filename_a = "test_a.bed"
        filename_b = "test_b.bed"
        intempdir() do
            write_intervals(filename_a, intervals_a)
            write_intervals(filename_b, intervals_b)
            @test check_intersection(filename_a, filename_b)
        end

    end

    @testset "GFF3 Parsing" begin
        get_bio_fmt_specimens()
        function check_gff3_parse(filename)
            # Reading from a stream
            for interval in GFF3Reader(open(filename))
            end

            # Reading from a regular file
            for interval in open(GFF3Reader, filename)
            end

            # in-place parsing
            stream = open(GFF3Reader, filename)
            entry = eltype(stream)()
            while !eof(stream)
                read!(stream, entry)
            end
            close(stream)

            return true
        end

        path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "GFF3")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            valid = get(specimen, "valid", true)
            if valid
                @test check_gff3_parse(joinpath(path, specimen["filename"]))
            else
                @test_throws Exception check_gff3_parse(joinpath(path, specimen["filename"]))
            end
        end
    end
end


@testset "BigBed" begin
    @testset "BED → BigBed → BED round-trip" begin
        path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "BED")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            if !get(specimen, "valid", true)
                continue
            end

            # BED → BigBed
            intervals = IntervalCollection(
                open(BEDReader, joinpath(path, specimen["filename"])))
            out = IOBuffer()
            write(BigBedWriter(out), intervals)

            # BigBed → BED
            seekstart(out)
            bb = BigBedReader(out)
            intervals2 = IntervalCollection(bb)

            @test intervals == intervals2
        end
    end

    @testset "BigBed Intersection" begin
        n = 10000
        srand(1234)
        chroms = ["one", "two", "three", "four", "five"]
        intervals = IntervalCollection(
            [Interval{BEDMetadata}(i.seqname, i.first, i.last, STRAND_NA, BEDMetadata())
             for i in random_intervals(chroms, 1000000, n)], true)

        # convert to bigbed in memory
        out = IOBuffer()
        write(BigBedWriter(out), intervals)
        seekstart(out)
        bb = BigBedReader(out)

        # intersection queries
        num_queries = 1000
        queries = random_intervals(chroms, 1000000, num_queries)

        @test all(Bool[IntervalCollection(collect(BEDInterval, intersect(intervals, query))) ==
                       IntervalCollection(collect(BEDInterval, intersect(bb, query)))
                       for query in queries])
    end

    # TODO: test summary information against output from kent's bigBedSummary
end

end # module TestIntervals
