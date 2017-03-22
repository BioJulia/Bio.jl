# Coverage
# ========
#
# Coverage of intervals.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
    coverage(intervals)

Compute the coverage of a collection of intervals and return an
`IntervalCollection` that contains run-length encoded coverage data.

For example, given intervals like:

    [------]     [------------]
       [---------------]

This function would return a new set of disjoint intervals with annotated
coverage like:

    [1][-2-][-1-][--2--][--1--]
"""
function coverage(stream, seqname_isless::Function=isless)
    cov = IntervalCollection{UInt32}()
    lasts = Int64[]

    stream_state = start(stream)
    if done(stream, stream_state)
        return cov
    end

    current_coverage = 0
    coverage_seqname = ""
    coverage_first = 0
    last_interval_first = 0
    interval, stream_state = next(stream, stream_state)
    while true
        if interval.seqname != coverage_seqname
            coverage_process_lasts_heap!(cov, current_coverage, coverage_seqname,
                                         coverage_first, lasts)
            if !(isempty(coverage_seqname) || seqname_isless(coverage_seqname, interval.seqname))
                error("Intervals must be sorted to compute coverage.")
            end

            coverage_seqname = interval.seqname
            current_coverage = 0
            coverage_first = 0
            last_interval_first = 0
        end

        if interval.first < last_interval_first
            error("Intervals must be sorted to compute coverage.")
        end

        if !isempty(lasts) && lasts[1] < first(interval)
            pos = Collections.heappop!(lasts)
            if first(interval) == pos + 1
                Collections.heappush!(lasts, last(interval))
                if done(stream, stream_state)
                    break
                end
                last_interval_first = first(interval)
                interval, stream_state = next(stream, stream_state)
            elseif pos == coverage_first - 1
                current_coverage -= 1
            else
                @assert pos >= coverage_first
                push!(cov, Interval{UInt32}(coverage_seqname, coverage_first,
                                            pos, STRAND_BOTH, current_coverage))
                current_coverage -= 1
                coverage_first = pos + 1
            end
        else
            if coverage_first == 0
                coverage_first = first(interval)
                current_coverage = 1
            elseif coverage_first == first(interval)
                current_coverage += 1
            else
                if current_coverage > 0
                    push!(cov, Interval{UInt32}(coverage_seqname, coverage_first,
                                                first(interval) - 1, STRAND_BOTH,
                                                current_coverage))
                end
                current_coverage += 1
                coverage_first = first(interval)
            end

            Collections.heappush!(lasts, last(interval))
            if done(stream, stream_state)
                break
            end
            last_interval_first = first(interval)
            interval, stream_state = next(stream, stream_state)
        end
    end

    coverage_process_lasts_heap!(cov, current_coverage, coverage_seqname,
                                 coverage_first, lasts)

    return cov
end

function coverage(ic::IntervalCollection)
    return coverage(ic, isless)
end

# Helper function for coverage. Process remaining interval end points after
# all intervals have been read.
function coverage_process_lasts_heap!(cov::IntervalCollection{UInt32},
                                      current_coverage, coverage_seqname,
                                      coverage_first, lasts)
    while !isempty(lasts)
        pos = Collections.heappop!(lasts)
        if pos == coverage_first - 1
            current_coverage -= 1
        else
            @assert pos >= coverage_first
            push!(cov, Interval{UInt32}(coverage_seqname, coverage_first,
                                        pos, STRAND_BOTH, current_coverage))
            current_coverage -= 1
            coverage_first = pos + 1
        end
    end
    @assert current_coverage == 0
end
