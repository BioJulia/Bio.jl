#
# Generic Feature Format Version 3 (GFF3)
# =======================================
#
# http://www.sequenceontology.org/gff3.shtml

immutable GFF3 <: FileFormat end

type GFF3Metadata
    source::StringField
    # use "kind" instead of "type" since "type" is a keyword of Julia
    kind::StringField
    score::Nullable{Float64}
    phase::Nullable{Int}
    attributes::StringField
end

function GFF3Metadata()
    return GFF3Metadata("", "", NaN, 0, "")
end

function Base.copy(metadata::GFF3Metadata)
    return GFF3Metadata(
        copy(metadata.source),
        copy(metadata.kind),
        metadata.score,
        metadata.phase,
        copy(metadata.attributes)
    )
end

"An `Interval` with associated metadata from a GFF3 file"
typealias GFF3Interval Interval{GFF3Metadata}

%%{
    machine _gff3parser;

    action finish_match {
        yield = true
        Ragel.@anchor!
        fbreak;
    }

    action count_line { input.state.linenum += 1 }
    action anchor { Ragel.anchor!(state, p) }

    action directive {
        directive = Ragel.@ascii_from_anchor!
        if startswith(directive, "##gff-version")
            # ##gff-version 3.2.1
            input.version = VersionNumber(split(directive, r"\s+")[2])
        elseif startswith(directive, "##sequence-region")
            # ##sequence-region seqid start end
            vals = split(directive, r"\s+")
            push!(input.sequence_regions,
                Interval(vals[2], parse(Int, vals[3]), parse(Int, vals[4])))
        else
            # TODO: record other directives
        end
    }

    # interval
    action seqname { Ragel.@copy_from_anchor!(output.seqname) }
    action start   { output.first = Ragel.@int64_from_anchor! }
    action end     { output.last = Ragel.@int64_from_anchor! }
    action strand  { output.strand = convert(Strand, (Ragel.@char)) }

    # metadata
    action source     { Ragel.@copy_from_anchor!(output.metadata.source) }
    action kind       { Ragel.@copy_from_anchor!(output.metadata.kind) }
    action score      { output.metadata.score = Nullable(Ragel.@float64_from_anchor!) }
    action nullscore  { output.metadata.score = Nullable{Float64}() }
    action phase      { output.metadata.phase = Nullable(Ragel.@int64_from_anchor!) }
    action nullphase  { output.metadata.phase = Nullable{Int}() }
    action attributes { Ragel.@copy_from_anchor!(output.metadata.attributes) }

    newline        = '\r'? '\n' >count_line;
    hspace         = [ \t\v];
    blankline      = hspace* newline;
    #floating_point = [\-+]?[0-9]* '.'? [0-9]+([eE][\-+]?[0-9]+)?;
    floating_point = [\-+] digit* '.'? digit+;

    comment   = '#' (any - newline)* newline;
    directive = ("##" (any - newline)* newline) >anchor %directive;

    seqname    = [ -~]* >anchor %seqname;
    source     = [ -~]* >anchor %source;
    kind       = [ -~]* >anchor %kind;
    start      = digit+ >anchor %start;
    end        = digit+ >anchor %end;
    score      = ((floating_point %score) | ('.' %nullscore)) >anchor;
    strand     = [+\-\.?] >strand;
    phase      = (([1-3] %phase) | ('.' %nullphase)) >anchor;
    attributes = (any - newline)* >anchor %attributes;

    gff3_entry = seqname '\t' source '\t' kind  '\t' start '\t' end '\t'
                 score   '\t' strand '\t' phase '\t' attributes newline
                 (comment | blankline)*;

    main := directive* comment* (gff3_entry %finish_match)*;
}%%

%% write data;


type GFF3Parser <: AbstractParser
    state::Ragel.State
    version::VersionNumber
    sequence_regions::Vector{Interval{Void}}

    function GFF3Parser(input::BufferedInputStream)
        %% write init;

        return new(Ragel.State(cs, input), VersionNumber(0), [])
    end
end

function Intervals.metadatatype(::GFF3Parser)
    return GFF3Metadata
end

function Base.eltype(::Type{GFF3Parser})
    return GFF3Interval
end

function open(input::BufferedInputStream, ::Type{GFF3})
    return GFF3Parser(input)
end

Ragel.@generate_read_fuction("_gff3parser", GFF3Parser, GFF3Interval,
    begin
        %% write exec;
    end)
