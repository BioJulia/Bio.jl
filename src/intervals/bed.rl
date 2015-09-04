

immutable BED <: FileFormat end


"""Metadata for BED interval records"""
type BEDMetadata
    used_fields::Int # how many of the first n fields are used
    name::StringField
    score::Int
    thick_first::Int
    thick_last::Int
    item_rgb::RGB{Float32}
    block_count::Int
    block_sizes::Vector{Int}
    block_firsts::Vector{Int}
end


function BEDMetadata()
    return BEDMetadata(0, StringField(), 0, 0, 0, RGB{Float32}(0.0, 0.0, 0.0),
                       0, Int[], Int[])
end


function Base.copy(metadata::BEDMetadata)
    return BEDMetadata(
        metadata.used_fields, copy(metadata.name),
        metadata.score, metadata.thick_first, metadata.thick_last,
        metadata.item_rgb, metadata.block_count,
        metadata.block_sizes[1:metadata.block_count],
        metadata.block_firsts[1:metadata.block_count])
end


function (==)(a::BEDMetadata, b::BEDMetadata)
    if a.used_fields != b.used_fields
        return false
    end

    n = a.used_fields
    ans = (n < 1 || a.name == b.name) &&
          (n < 2 || a.score == b.score) &&
          (n < 3 || a.thick_first == b.thick_first) &&
          (n < 4 || a.thick_last == b.thick_first) &&
          (n < 5 || a.item_rgb == b.item_rgb) &&
          (n < 6 || a.block_count == b.block_count)
    if !ans
        return false
    end

    if n >= 7
        for i in 1:a.block_count
            if a.block_sizes[i] != b.block_sizes[i]
                return false
            end
        end
    end

    if n >= 8
        for i in 1:a.block_count
            if a.block_sizes[i] != b.block_sizes[i]
                return false
            end
        end
    end

    return true
end

# TODO
#function show(io::IO, metadata::BEDMetadata)
#end


"An `Interval` with associated metadata from a BED file"
typealias BEDInterval Interval{BEDMetadata}


module BEDParserImpl

import Bio.Ragel, Bio.Intervals
using Bio: AbstractParser, StringField
using Bio.Intervals: Strand, STRAND_NA, BED, BEDInterval, BEDMetadata
using BufferedStreams, Switch, Compat, Colors


%%{
    machine bed;

    action finish_match {
        input.block_size_idx = 1
        input.block_first_idx = 1
        output.metadata.used_fields = 0

        yield = true
        # // fbreak causes will cause the pushmark action for the next seqname
        # // to be skipped, so we do it here TODO: Is this still the case????
        Ragel.anchor!(state, p)
        fbreak;
    }

    action count_line { input.state.linenum += 1 }
    action anchor { Ragel.anchor!(state, p) }
    action optional_field { output.metadata.used_fields += 1 }

    action seqname     { Ragel.@copy_from_anchor!(output.seqname) }
    action first       { output.first = 1 + Ragel.@int64_from_anchor! }
    action last        { output.last = Ragel.@int64_from_anchor! }
    action name        { Ragel.@copy_from_anchor!(output.metadata.name) }
    action score       { output.metadata.score = Ragel.@int64_from_anchor! }
    action strand      { output.strand = convert(Strand, (Ragel.@char)) }
    action thick_first { output.metadata.thick_first = 1 + Ragel.@int64_from_anchor! }
    action thick_last  { output.metadata.thick_last = Ragel.@int64_from_anchor!  }
    action item_rgb_r  { input.red = input.green = input.blue = (Ragel.@int64_from_anchor!) / 255.0 }
    action item_rgb_g  { input.green = (Ragel.@int64_from_anchor!) / 255.0 }
    action item_rgb_b  { input.blue = (Ragel.@int64_from_anchor!) / 255.0 }
    action item_rgb    { output.metadata.item_rgb = RGB{Float32}(input.red, input.green, input.blue) }
    action block_count {
        output.metadata.block_count = Ragel.@int64_from_anchor!

        if (output.metadata.block_count > length(output.metadata.block_sizes))
            resize!(output.metadata.block_sizes, output.metadata.block_count)
        end

        if (output.metadata.block_count > length(output.metadata.block_firsts))
            resize!(output.metadata.block_firsts, output.metadata.block_count)
        end
    }

    action block_size {
        if input.block_size_idx > length(output.metadata.block_sizes)
            error("More size blocks encountered than BED block count field suggested.")
        end
        output.metadata.block_sizes[input.block_size_idx] = Ragel.@int64_from_anchor!
        input.block_size_idx += 1
    }

    action block_first {
        if input.block_first_idx > length(output.metadata.block_firsts)
            error("More start blocks encountered than BED block count field suggested.")
        end
        output.metadata.block_firsts[input.block_first_idx] = Ragel.@int64_from_anchor!
        input.block_first_idx += 1
    }

    newline      = '\r'? '\n'     >count_line;
    hspace       = [ \t\v];
    blankline    = hspace* newline;

    seqname      = [ -~]*   >anchor     %seqname;
    first        = digit+   >anchor     %first;
    last         = digit+   >anchor     %last;
    name         = [ -~]*   >anchor     %name  %optional_field;
    score        = digit+   >anchor     %score %optional_field;
    strand       = [+\-\.?] >strand %optional_field;
    thick_first  = digit+   >anchor     %thick_first %optional_field;
    thick_last   = digit+   >anchor     %thick_last  %optional_field;

    item_rgb_r   = digit+   >anchor     %item_rgb_r;
    item_rgb_g   = digit+   >anchor     %item_rgb_g;
    item_rgb_b   = digit+   >anchor     %item_rgb_b;
    item_rgb     = item_rgb_r (hspace* ',' hspace* item_rgb_g hspace* ',' hspace* item_rgb_b)? %item_rgb %optional_field;

    block_count  = digit+   >anchor    %block_count %optional_field;

    block_size   = digit+   >anchor    %block_size;
    block_sizes  = block_size (',' block_size)* ','? %optional_field;

    block_first  = digit+   >anchor    %block_first;
    block_firsts = block_first (',' block_first)* ','? %optional_field;

    bed_entry = seqname '\t' first '\t' last (
                    '\t' name ( '\t' score ( '\t' strand ( '\t' thick_first (
                    '\t' thick_last ( '\t' item_rgb ( '\t' block_count (
                    '\t' block_sizes ( '\t' block_firsts )? )? )? )? )? )? )? )?
                    )?
                newline blankline*;

    main := blankline* (bed_entry %finish_match)*;
}%%


%% write data;


type BEDParser <: AbstractParser
    state::Ragel.State

    # intermediate values used during parsing
    red::Float32
    green::Float32
    blue::Float32
    block_size_idx::Int
    block_first_idx::Int

    function BEDParser(input::BufferedInputStream)
        %% write init;

        return new(Ragel.State(cs, input), 0.0, 0.0, 0.0, 1, 1)
    end
end


function Intervals.metadatatype(::BEDParser)
    return BEDMetadata
end


function Base.eltype(::Type{BEDParser})
    return BEDInterval
end


function Base.open(input::BufferedInputStream, ::Type{BED})
    return BEDParser(input)
end


Ragel.@generate_read_fuction("bed", BEDParser, BEDInterval,
    begin
        %% write exec;
    end)



end # module BEDParserImpl


# TODO: Rewrite this stuff

"""
Write a BEDInterval in BED format.
"""
function write(out::IO, interval::BEDInterval)
    print(out, interval.seqname, '\t', interval.first - 1, '\t', interval.last)
    write_optional_fields(out, interval)
    write(out, '\n')
end


function write_optional_fields(out::IO, interval::BEDInterval, leadingtab::Bool=true)
    if !isnull(interval.metadata.name)
        if leadingtab
            write(out, '\t')
        end
        write(out, get(interval.metadata.name))
    else return end

    if !isnull(interval.metadata.score)
        print(out, '\t', get(interval.metadata.score))
    else return end

    if interval.strand != STRAND_NA
        print(out, '\t', interval.strand)
    else return end

    if !isnull(interval.metadata.thick_first)
        print(out, '\t', get(interval.metadata.thick_first) - 1)
    else return end

    if !isnull(interval.metadata.thick_last)
        print(out, '\t', get(interval.metadata.thick_last))
    else return end

    if !isnull(interval.metadata.item_rgb)
        item_rgb = get(interval.metadata.item_rgb)
        print(out, '\t',
              round(Int, 255 * item_rgb.r), ',',
              round(Int, 255 * item_rgb.g), ',',
              round(Int, 255 * item_rgb.b))
    else return end

    if !isnull(interval.metadata.block_count)
        print(out, '\t', get(interval.metadata.block_count))
    else return end

    if !isnull(interval.metadata.block_sizes)
        block_sizes = get(interval.metadata.block_sizes)
        if !isempty(block_sizes)
            print(out, '\t', block_sizes[1])
            for i in 2:length(block_sizes)
                print(out, ',', block_sizes[i])
            end
        end
    else return end

    if !isnull(interval.metadata.block_firsts)
        block_firsts = get(interval.metadata.block_firsts)
        if !isempty(block_firsts)
            print(out, '\t', block_firsts[1] - 1)
            for i in 2:length(block_firsts)
                print(out, ',', block_firsts[i] - 1)
            end
        end
    end
end


function IntervalCollection(interval_stream::BEDParserImpl.BEDParser)
    intervals = collect(BEDInterval, interval_stream)
    return IntervalCollection{BEDMetadata}(intervals, true)
end


