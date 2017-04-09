# SAM Reader
# =========

type Reader <: Bio.IO.AbstractReader
    state::Bio.Ragel.State
    header::Header

    function Reader(input::BufferedStreams.BufferedInputStream)
        reader = new(Bio.Ragel.State(sam_header_machine.start_state, input), Header())
        readheader!(reader)
        reader.state.cs = sam_body_machine.start_state
        return reader
    end
end

"""
    SAM.Reader(input::IO)

Create a data reader of the SAM file format.

# Arguments
* `input`: data source
"""
function Reader(input::IO)
    return Reader(BufferedStreams.BufferedInputStream(input))
end

function Bio.IO.stream(reader::Reader)
    return reader.state.stream
end

"""
    header(reader::Reader)::Header

Get the header of `reader`.
"""
function header(reader::Reader)::Header
    return reader.header
end

function Bio.header(reader::Reader)
    return header(reader)
end

function Base.eltype(::Type{Reader})
    return Record
end

# file   = header . body
# header = metainfo*
# body   = record*
info("compiling SAM")
const sam_metainfo_machine, sam_record_machine, sam_header_machine, sam_body_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    any = Automa.RegExp.any

    metainfo = let
        tag = re"[A-Z][A-Z]" \ cat("CO")
        tag.actions[:enter] = [:mark1]
        tag.actions[:exit]  = [:metainfo_tag]

        dict = let
            key = re"[A-Za-z][A-Za-z0-9]"
            key.actions[:enter] = [:mark2]
            key.actions[:exit]  = [:metainfo_dict_key]
            val = re"[ -~]+"
            val.actions[:enter] = [:mark2]
            val.actions[:exit]  = [:metainfo_dict_val]
            keyval = cat(key, ':', val)

            cat(keyval, rep(cat('\t', keyval)))
        end
        dict.actions[:enter] = [:mark1]
        dict.actions[:exit]  = [:metainfo_val]

        co = cat("CO")
        co.actions[:enter] = [:mark1]
        co.actions[:exit]  = [:metainfo_tag]

        comment = re"[^\r\n]*"
        comment.actions[:enter] = [:mark1]
        comment.actions[:exit]  = [:metainfo_val]

        cat('@', alt(cat(tag, '\t', dict), cat(co, '\t', comment)))
    end
    metainfo.actions[:enter] = [:anchor]
    metainfo.actions[:exit]  = [:metainfo]

    record = let
        qname = re"[!-?A-~]+"
        qname.actions[:enter] = [:mark]
        qname.actions[:exit]  = [:record_qname]

        flag = re"[0-9]+"
        flag.actions[:enter] = [:mark]
        flag.actions[:exit]  = [:record_flag]

        rname = re"\*|[!-()+-<>-~][!-~]*"
        rname.actions[:enter] = [:mark]
        rname.actions[:exit]  = [:record_rname]

        pos = re"[0-9]+"
        pos.actions[:enter] = [:mark]
        pos.actions[:exit]  = [:record_pos]

        mapq = re"[0-9]+"
        mapq.actions[:enter] = [:mark]
        mapq.actions[:exit]  = [:record_mapq]

        cigar = re"\*|([0-9]+[MIDNSHPX=])+"
        cigar.actions[:enter] = [:mark]
        cigar.actions[:exit]  = [:record_cigar]

        rnext = re"\*|=|[!-()+-<>-~][!-~]*"
        rnext.actions[:enter] = [:mark]
        rnext.actions[:exit]  = [:record_rnext]

        pnext = re"[0-9]+"
        pnext.actions[:enter] = [:mark]
        pnext.actions[:exit]  = [:record_pnext]

        tlen = re"[-+]?[0-9]+"
        tlen.actions[:enter] = [:mark]
        tlen.actions[:exit]  = [:record_tlen]

        seq = re"\*|[A-Za-z=.]+"
        seq.actions[:enter] = [:mark]
        seq.actions[:exit]  = [:record_seq]

        qual = re"[!-~]+"
        qual.actions[:enter] = [:mark]
        qual.actions[:exit]  = [:record_qual]

        field = let
            tag = re"[A-Za-z][A-Za-z0-9]"
            val = alt(
                re"A:[!-~]",
                re"i:[-+]?[0-9]+",
                re"f:[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?",
                re"Z:[ !-~]*",
                re"H:([0-9A-F][0-9A-F])*",
                re"B:[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+")

            cat(tag, ':', val)
        end
        field.actions[:enter] = [:mark]
        field.actions[:exit]  = [:record_field]

        cat(
            qname, '\t',
            flag,  '\t',
            rname, '\t',
            pos,   '\t',
            mapq,  '\t',
            cigar, '\t',
            rnext, '\t',
            pnext, '\t',
            tlen,  '\t',
            seq,   '\t',
            qual,
            rep(cat('\t', field)))
    end
    record.actions[:enter] = [:anchor]
    record.actions[:exit]  = [:record]

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(re"\r?", lf)
    end

    header′ = rep(cat(metainfo, newline))
    header′.actions[:exit] = [:header]
    header = cat(header′, opt(any() \ cat('@')))  # look ahead

    body = rep(cat(record, newline))

    return map(Automa.compile, (metainfo, record, header, body))
end)()

const sam_metainfo_actions = Dict(
    :metainfo_tag => :(record.tag = (mark1:p-1) - offset),
    :metainfo_val => :(record.val = (mark1:p-1) - offset),
    :metainfo_dict_key => :(push!(record.dictkey, (mark2:p-1) - offset)),
    :metainfo_dict_val => :(push!(record.dictval, (mark2:p-1) - offset)),
    :metainfo => quote
        Bio.ReaderHelper.resize_and_copy!(record.data, data, offset+1:p-1)
        record.filled = (offset+1:p-1) - offset
    end,
    :anchor => :(),
    :mark1  => :(mark1 = p),
    :mark2  => :(mark2 = p))
eval(
    Bio.ReaderHelper.generate_index_function(
        MetaInfo,
        sam_metainfo_machine,
        :(mark1 = mark2 = offset = 0),
        sam_metainfo_actions))
eval(
    Bio.ReaderHelper.generate_readheader_function(
        Reader,
        MetaInfo,
        sam_header_machine,
        :(mark1 = mark2 = offset = 0),
        merge(sam_metainfo_actions, Dict(
            :metainfo => quote
                Bio.ReaderHelper.resize_and_copy!(record.data, data, Bio.ReaderHelper.upanchor!(stream):p-1)
                record.filled = (offset+1:p-1) - offset
                @assert isfilled(record)
                push!(reader.header.metainfo, record)
                Bio.ReaderHelper.ensure_margin!(stream)
                record = MetaInfo()
            end,
            :header => :(finish_header = true; @escape),
            :countline => :(linenum += 1),
            :anchor => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1))),
        quote
            if !eof(stream)
                stream.position -= 1  # cancel look-ahead
            end
        end))

const sam_record_actions = Dict(
    :record_qname => :(record.qname = (mark:p-1) - offset),
    :record_flag  => :(record.flag  = (mark:p-1) - offset),
    :record_rname => :(record.rname = (mark:p-1) - offset),
    :record_pos   => :(record.pos   = (mark:p-1) - offset),
    :record_mapq  => :(record.mapq  = (mark:p-1) - offset),
    :record_cigar => :(record.cigar = (mark:p-1) - offset),
    :record_rnext => :(record.rnext = (mark:p-1) - offset),
    :record_pnext => :(record.pnext = (mark:p-1) - offset),
    :record_tlen  => :(record.tlen  = (mark:p-1) - offset),
    :record_seq   => :(record.seq   = (mark:p-1) - offset),
    :record_qual  => :(record.qual  = (mark:p-1) - offset),
    :record_field => :(push!(record.fields, (mark:p-1) - offset)),
    :record       => quote
        Bio.ReaderHelper.resize_and_copy!(record.data, data, 1:p-1)
        record.filled = (offset+1:p-1) - offset
    end,
    :anchor       => :(),
    :mark         => :(mark = p))
eval(
    Bio.ReaderHelper.generate_index_function(
        Record,
        sam_record_machine,
        :(mark = offset = 0),
        sam_record_actions))
eval(
    Bio.ReaderHelper.generate_read_function(
        Reader,
        sam_body_machine,
        :(mark = offset = 0),
        merge(sam_record_actions, Dict(
            :record    => quote
                Bio.ReaderHelper.resize_and_copy!(record.data, data, Bio.ReaderHelper.upanchor!(stream):p-1)
                record.filled = (offset+1:p-1) - offset
                found_record = true
                @escape
            end,
            :countline => :(linenum += 1),
            :anchor    => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))))
