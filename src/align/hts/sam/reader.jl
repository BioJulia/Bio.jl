# SAM Reader
# ==========

"""
    SAMReader(input::IO)

Create a data reader of the SAM file format.

# Arguments
* `input`: data source
"""
type SAMReader <: Bio.IO.AbstractReader
    state::Bio.Ragel.State
    header::SAMHeader

    function SAMReader(input::BufferedStreams.BufferedInputStream)
        reader = new(Bio.Ragel.State(sam_header_machine.start_state, input), SAMHeader())
        readheader!(reader)
        reader.state.cs = sam_body_machine.start_state
        return reader
    end
end

function SAMReader(input::IO)
    return SAMReader(BufferedStreams.BufferedInputStream(input))
end

function Bio.IO.stream(reader::SAMReader)
    return reader.state.stream
end

function header(reader::SAMReader)
    return reader.header
end

function Base.eltype(::Type{SAMReader})
    return SAMRecord
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
    :metainfo => :(),
    :anchor => :(),
    :mark1  => :(mark1 = p),
    :mark2  => :(mark2 = p))

eval(Bio.ReaderHelper.generate_index_function(SAMMetaInfo, sam_metainfo_machine, sam_metainfo_actions))

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
    :record       => :(),
    :anchor       => :(),
    :mark         => :(mark = p))

eval(Bio.ReaderHelper.generate_index_function(SAMRecord, sam_record_machine, sam_record_actions))

const sam_header_actions = merge(sam_metainfo_actions, Dict(
    :metainfo => quote
        record.data = data[Bio.ReaderHelper.upanchor!(stream):p-1]
        record.filled = true
        push!(reader.header.metainfo, record)
        record = SAMMetaInfo()
    end,
    :header => :(finish_header = true; @escape),
    :countline => :(linenum += 1),
    :anchor => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))

function readheader!(reader::SAMReader)
    _readheader!(reader, reader.state)
end

@eval function _readheader!(reader::SAMReader, state::Bio.Ragel.State)
    stream = state.stream
    Bio.ReaderHelper.ensure_margin!(stream)
    cs = state.cs
    linenum = state.linenum
    data = stream.buffer
    p = stream.position
    p_end = stream.available
    p_eof = -1
    offset = mark1 = mark2 = 0
    finish_header = false
    record = SAMMetaInfo()
 
    while true
        $(Automa.generate_exec_code(sam_header_machine, actions=sam_header_actions, code=:table))
 
        state.cs = cs
        state.finished = cs == 0
        state.linenum = linenum
        stream.position = p

        if cs < 0
            error("SAM file format error on line ", linenum)
        elseif finish_header
            #upanchor!(stream)
            if !eof(stream)
                stream.position -= 1  # cancel look ahead
            end
            break
        #elseif cs == 0
        #    throw(EOFError())
        elseif p > p_eof ≥ 0
            error("incomplete SAM input on line ", linenum)
        else
            hits_eof = BufferedStreams.fillbuffer!(stream) == 0
            p = stream.position
            p_end = stream.available
            if hits_eof
                p_eof = p_end
            end
        end
    end
end

const sam_body_actions = merge(sam_record_actions, Dict(
    :record    => :(found_record = true; @escape),
    :countline => :(linenum += 1),
    :anchor    => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))

eval(Bio.ReaderHelper.generate_read_function(SAMReader, sam_body_machine, sam_body_actions))
