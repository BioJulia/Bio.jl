# GFF3 Reader
# ===========

"""
    GFF3Reader(input::IO; save_directives::Bool=false)

Create a reader for data in GFF3 format.

# Arguments:
* `input`: data source
* `save_directives=false`: if true, store directive lines, which can be accessed
  with the `directives` function
"""
type Reader <: Bio.IO.AbstractReader
    state::Bio.Ragel.State
    save_directives::Bool
    targets::Vector{Symbol}
    found_fasta::Bool
    directives::Vector{Record}
    directive_count::Int
    preceding_directive_count::Int

    function Reader(input::BufferedStreams.BufferedInputStream,
                    save_directives::Bool=false,
                    skip_features::Bool=false, skip_directives::Bool=true, skip_comments::Bool=true)
        targets = Symbol[]
        if !skip_features
            push!(targets, :feature)
        end
        if !skip_directives
            push!(targets, :directive)
        end
        if !skip_comments
            push!(targets, :comment)
        end
        return new(Bio.Ragel.State(body_machine.start_state, input), save_directives, targets, false, Record[], 0, 0)
    end
end

function Reader(input::IO;
                save_directives::Bool=false,
                skip_features::Bool=false, skip_directives::Bool=true, skip_comments::Bool=true)
    return Reader(BufferedStreams.BufferedInputStream(input), save_directives, skip_features, skip_directives, skip_comments)
end

function Base.eltype(::Type{Reader})
    return Record
end

function Bio.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.eof(reader::Reader)
    return reader.state.finished || eof(reader.state.stream)
end

function Base.close(reader::Reader)
    # make trailing directives accessable
    reader.directive_count = reader.preceding_directive_count
    reader.preceding_directive_count = 0
    close(Bio.IO.stream(reader))
end

function IntervalCollection(reader::Reader)
    intervals = collect(Interval{Record}, reader)
    return IntervalCollection(intervals, true)
end


"""
Return all directives that preceded the last GFF entry parsed as an array of
strings.

Directives at the end of the file can be accessed by calling `close(reader)`
and then `directives(reader)`.
"""
function directives(reader::Reader)
    ret = String[]
    for i in endof(reader.directives)-reader.directive_count+1:endof(reader.directives)
        push!(ret, content(reader.directives[i]))
    end
    return ret
end

"""
Return true if the GFF3 stream is at its end and there is trailing FASTA data.
"""
function hasfasta(reader::Reader)
    if eof(reader)
        return reader.found_fasta
    else
        error("GFF3 file must be read until the end before any FASTA sequences can be accessed")
    end
end

"""
Return a Bio.Seq.FASTA.Reader initialized to parse trailing FASTA data.

Throws an exception if there is no trailing FASTA, which can be checked using
`hasfasta`.
"""
function getfasta(reader::Reader)
    if !hasfasta(reader)
        error("GFF3 file has no FASTA data")
    end
    return Bio.Seq.FASTA.Reader(reader.state.stream)
end

info("compiling GFF3")
const record_machine, body_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    rep1 = Automa.RegExp.rep1
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt

    feature = let
        seqid = re"[a-zA-Z0-9.:^*$@!+_?\-|%]*"
        seqid.actions[:enter] = [:mark]
        seqid.actions[:exit]  = [:feature_seqid]

        source = re"[ -~]*"
        source.actions[:enter] = [:mark]
        source.actions[:exit]  = [:feature_source]

        type_ = re"[ -~]*"
        type_.actions[:enter] = [:mark]
        type_.actions[:exit]  = [:feature_type_]

        start = re"[0-9]+|\."
        start.actions[:enter] = [:mark]
        start.actions[:exit]  = [:feature_start]

        end_ = re"[0-9]+|\."
        end_.actions[:enter] = [:mark]
        end_.actions[:exit]  = [:feature_end_]

        score = re"[ -~]*[0-9][ -~]*|\."
        score.actions[:enter] = [:mark]
        score.actions[:exit]  = [:feature_score]

        strand = re"[+\-?]|\."
        strand.actions[:enter] = [:feature_strand]

        phase = re"[0-2]|\."
        phase.actions[:enter] = [:feature_phase]

        attributes = let
            char = re"[^=;,\t\r\n]"
            key = rep1(char)
            key.actions[:enter] = [:mark]
            key.actions[:exit]  = [:feature_attribute_key]
            val = rep(char)
            attr = cat(key, '=', val, rep(cat(',', val)))

            cat(rep(cat(attr, ';')), opt(attr))
        end

        cat(seqid,  '\t',
            source, '\t',
            type_,  '\t',
            start,  '\t',
            end_,   '\t',
            score,  '\t',
            strand, '\t',
            phase,  '\t',
            attributes)
    end
    feature.actions[:exit] = [:feature]

    directive = re"##[^\r\n]*"
    directive.actions[:exit] = [:directive]

    comment = re"#([^#\r\n][^\r\n]*)?"
    comment.actions[:exit] = [:comment]

    record = alt(feature, directive, comment)
    record.actions[:enter] = [:anchor]
    record.actions[:exit]  = [:record]

    blank = re"[ \t]*"

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(opt('\r'), lf)
    end

    body = rep(cat(alt(record, blank), newline))
    body.actions[:exit] = [:body]

    # look-ahead of the beginning of FASTA
    body′ = cat(body, opt('>'))

    return map(Automa.compile, (record, body′))
end)()

const record_actions = Dict(
    :feature_seqid   => :(record.seqid  = (mark:p-1) - offset),
    :feature_source  => :(record.source = (mark:p-1) - offset),
    :feature_type_   => :(record.type_  = (mark:p-1) - offset),
    :feature_start   => :(record.start  = (mark:p-1) - offset),
    :feature_end_    => :(record.end_   = (mark:p-1) - offset),
    :feature_score   => :(record.score  = (mark:p-1) - offset),
    :feature_strand  => :(record.strand = p - offset),
    :feature_phase   => :(record.phase  = p - offset),
    :feature_attribute_key => :(push!(record.attribute_keys, (mark:p-1) - offset)),
    :feature         => :(record.kind = :feature),
    :directive       => :(record.kind = :directive),
    :comment         => :(record.kind = :comment),
    :record          => quote
        Bio.ReaderHelper.resize_and_copy!(record.data, data, 1:p-1)
        record.filled = (offset+1:p-1) - offset
    end,
    :anchor          => :(),
    :mark            => :(mark = p))
eval(
    Bio.ReaderHelper.generate_index_function(
        Record,
        record_machine,
        quote
            mark = offset = 0
        end,
        record_actions))
eval(
    Bio.ReaderHelper.generate_read_function(
        Reader,
        body_machine,
        quote
            mark = offset = 0
        end,
        merge(record_actions, Dict(
            :record => quote
                Bio.ReaderHelper.resize_and_copy!(record.data, data, Bio.ReaderHelper.upanchor!(stream):p-1)
                record.filled = (offset+1:p-1) - offset
                if isfeature(record)
                    reader.directive_count = reader.preceding_directive_count
                    reader.preceding_directive_count = 0
                elseif isdirective(record)
                    reader.preceding_directive_count += 1
                    if reader.save_directives
                        push!(reader.directives, copy(record))
                    end
                    if is_fasta_directive(record)
                        reader.found_fasta = true
                        reader.state.finished = true
                    end
                end
                if record.kind ∈ reader.targets
                    found_record = true
                    @escape
                end
            end,
            :body => quote
                if data[p] == UInt8('>')
                    reader.found_fasta = true
                    reader.state.finished = true
                    # HACK: any better way?
                    cs = 0
                    @goto exit
                end
            end,
            :countline => :(linenum += 1),
            :anchor    => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))))
