# GFF3 Reader
# ===========

type Reader <: Bio.IO.AbstractReader
    state::Bio.Ragel.State
    save_directives::Bool
    found_fasta::Bool
    directives::Vector{Record}

    function Reader(input::BufferedStreams.BufferedInputStream, save_directives::Bool=false)
        return new(Bio.Ragel.State(body_machine.start_state, input), save_directives, false, Record[])
    end
end

function Reader(input::IO; save_directives::Bool=false)
    return Reader(BufferedStreams.BufferedInputStream(input), save_directives)
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
Return a FASTAReader initialized to parse trailing FASTA data.

Throws an exception if there is no trailing FASTA, which can be checked using
`hasfasta`.
"""
function getfasta(reader::Reader)
    if !hasfasta(reader)
        error("GFF3 file has no FASTA data")
    end
    return Bio.Seq.FASTAReader(reader.state.stream)
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

        typ = re"[ -~]*"
        typ.actions[:enter] = [:mark]
        typ.actions[:exit]  = [:feature_typ]

        start = re"[0-9]+|\."
        start.actions[:enter] = [:mark]
        start.actions[:exit]  = [:feature_start]

        stop = re"[0-9]+|\."
        stop.actions[:enter] = [:mark]
        stop.actions[:exit]  = [:feature_stop]

        score = re"[ -~]*[0-9][ -~]*|\."
        score.actions[:enter] = [:mark]
        score.actions[:exit]  = [:feature_score]

        strand = re"[+\-?]|\."
        strand.actions[:enter] = [:feature_strand]

        phase = re"[0-2]|\."
        phase.actions[:enter] = [:feature_phase]

        attributes = let
            char = re"[ -~]" \ re"[=;,]"
            key = rep1(char)
            key.actions[:enter] = [:mark]
            key.actions[:exit]  = [:feature_attribute_key]
            val = rep(char)
            attr = cat(key, '=', val, rep(cat(',', val)))

            opt(cat(attr, rep(cat(';', attr))))
        end

        cat(seqid,  '\t',
            source, '\t',
            typ,    '\t',
            start,  '\t',
            stop,   '\t',
            score,  '\t',
            strand, '\t',
            phase,  '\t',
            attributes)
    end
    feature.actions[:exit] = [:feature]

    directive = cat("##", re"[ -~]*")
    directive.actions[:exit] = [:directive]

    comment = cat('#', opt(cat(re"[ -~]" \ cat('#'), re"[ -~]*")))
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

    return map(Automa.compile, (record, body))
end)()

const record_actions = Dict(
    :feature_seqid   => :(record.seqid  = (mark:p-1) - offset),
    :feature_source  => :(record.source = (mark:p-1) - offset),
    :feature_typ     => :(record.typ    = (mark:p-1) - offset),
    :feature_start   => :(record.start  = (mark:p-1) - offset),
    :feature_stop    => :(record.stop   = (mark:p-1) - offset),
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
        record_actions))
eval(
    Bio.ReaderHelper.generate_read_function(
        Reader,
        body_machine,
        merge(record_actions, Dict(
            :record    => quote
                Bio.ReaderHelper.resize_and_copy!(record.data, data, Bio.ReaderHelper.upanchor!(stream):p-1)
                record.filled = (offset+1:p-1) - offset
                found_record = true
                if isdirective(record) && reader.save_directives
                    push!(reader.directive, copy(record))
                end
                if is_fasta_directive(record)
                    reader.found_fasta = true
                    reader.state.finished = true
                end
                @escape
            end,
            :countline => :(linenum += 1),
            :anchor    => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))))
