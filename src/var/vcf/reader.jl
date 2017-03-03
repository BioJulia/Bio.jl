# VCF Reader
# ==========

type VCFReader <: Bio.IO.AbstractReader
    state::Bio.Ragel.State
    header::VCFHeader

    function VCFReader(input::BufferedInputStream)
        reader = new(Bio.Ragel.State(vcf_header_machine.start_state, input), VCFHeader())
        readheader!(reader)
        reader.state.cs = vcf_body_machine.start_state
        return reader
    end
end

function VCFReader(input::IO)
    return VCFReader(BufferedInputStream(input))
end

function Base.eltype(::Type{VCFReader})
    return VCFRecord
end

function Bio.IO.stream(reader::VCFReader)
    return reader.state.stream
end

function header(reader::VCFReader)
    return reader.header
end

# VCF v4.3
info("compiling VCF")
const vcf_metainfo_machine, vcf_record_machine, vcf_header_machine, vcf_body_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    delim(x, sep) = opt(cat(x, rep(cat(sep, x))))

    # The 'fileformat' field is required and must be the first line.
    fileformat = let
        key = cat("fileformat")
        key.actions[:enter] = [:mark]
        key.actions[:exit]  = [:metainfo_key]

        version = re"[!-~]+"
        version.actions[:enter] = [:mark2]
        version.actions[:exit]  = [:metainfo_val]

        cat("##", key, '=', version)
    end
    fileformat.actions[:enter] = [:anchor]
    fileformat.actions[:exit]  = [:metainfo]

    # All kinds of meta-information line after 'fileformat' are handled here.
    metainfo = let
        key = re"[0-9A-Za-z_]+"
        key.actions[:enter] = [:mark]
        key.actions[:exit]  = [:metainfo_key]

        str = re"[ -;=-~][ -~]*"  # does not starts with '<'
        str.actions[:enter] = [:mark2]
        str.actions[:exit]  = [:metainfo_val]

        dict = let
            dictkey = re"[0-9A-Za-z_]+"
            dictkey.actions[:enter] = [:mark]
            dictkey.actions[:exit]  = [:metainfo_dict_key]

            dictval = let
                quoted   = cat('"', rep(alt(re"[ !#-[\]-~]", "\\\"", "\\\\")), '"')
                unquoted = rep(re"[ -~]" \ re"[\",>]")
                alt(quoted, unquoted)
            end
            dictval.actions[:enter] = [:mark]
            dictval.actions[:exit]  = [:metainfo_dict_val]

            cat('<', delim(cat(dictkey, '=', dictval), ','), '>')
        end
        dict.actions[:enter] = [:mark2]
        dict.actions[:exit]  = [:metainfo_val]

        cat("##", key, '=', alt(str, dict))
    end
    metainfo.actions[:enter] = [:anchor]
    metainfo.actions[:exit]  = [:metainfo]

    # The header line.
    header = let
        sampleID = re"[ -~]+"
        sampleID.actions[:enter] = [:mark]
        sampleID.actions[:exit]  = [:header_sampleID]

        cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", opt(re"\tFORMAT" * rep(re"\t" * sampleID)))
    end
    header.actions[:enter] = [:anchor]

    # Data lines (fixed fields and variable genotype fields).
    record = let
        chrom = re"[!-9;-~]+"  # no colon
        chrom.actions[:enter] = [:mark]
        chrom.actions[:exit]  = [:record_chrom]

        pos = re"[0-9]+|\."
        pos.actions[:enter] = [:mark]
        pos.actions[:exit]  = [:record_pos]

        id = let
            elm = re"[!-:<-~]+" \ cat('.')
            elm.actions[:enter] = [:mark]
            elm.actions[:exit]  = [:record_id]

            alt(delim(elm, ';'), '.')
        end

        ref = re"[!-~]+"
        ref.actions[:enter] = [:mark]
        ref.actions[:exit]  = [:record_ref]

        alt′ = let
            elm = re"[!-+--~]+" \ cat('.')
            elm.actions[:enter] = [:mark]
            elm.actions[:exit]  = [:record_alt]

            alt(delim(elm, ','), '.')
        end

        qual = re"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?|NaN|[-+]Inf|\."
        qual.actions[:enter] = [:mark]
        qual.actions[:exit]  = [:record_qual]

        filter = let
            elm = re"[!-:<-~]+" \ cat('.')
            elm.actions[:enter] = [:mark]
            elm.actions[:exit]  = [:record_filter]

            alt(delim(elm, ';'), '.')
        end

        info = let
            key = re"[A-Za-z_][0-9A-Za-z_.]*"
            key.actions[:enter] = [:mark]
            key.actions[:exit]  = [:record_info_key]

            val = opt(cat('=', re"[ -:<-~]+"))

            alt(delim(cat(key, val), ';'), '.')
        end

        format = let
            elm = re"[A-Za-z_][0-9A-Za-z_.]*"
            elm.actions[:enter] = [:mark]
            elm.actions[:exit]  = [:record_format]

            alt(delim(elm, ':'), '.')
        end

        genotype = let
            elm = re"[ -9;-~]+"  # no colon
            elm.actions[:enter] = [:mark]
            elm.actions[:exit]  = [:record_genotype_elm]

            delim(elm, ':')
        end
        genotype.actions[:enter] = [:record_genotype]

        cat(
            chrom,  '\t',
            pos,    '\t',
            id,     '\t',
            ref,    '\t',
            alt′,   '\t',
            qual,   '\t',
            filter, '\t',
            info,
            opt(cat('\t', format, rep(cat('\t', genotype)))))
    end
    record.actions[:enter] = [:anchor]
    record.actions[:exit]  = [:record]

    # A newline can be either a CR+LF or a LF.
    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(opt('\r'), lf)
    end

    # The VCF file format (header and body part).
    vcfheader = cat(
        fileformat, newline,
        rep(cat(metainfo, newline)),
        header, newline)
    vcfheader.actions[:final] = [:vcfheader]

    vcfbody = rep(cat(record, newline))

    return map(Automa.compile, (metainfo, record, vcfheader, vcfbody))
end)()

#= Debug
write("vcf_header.dot", Automa.dfa2dot(vcf_header_machine.dfa))
run(`dot -Tsvg -o vcf_header.svg vcf_header.dot`)
write("vcf_body.dot", Automa.dfa2dot(vcf_body_machine.dfa))
run(`dot -Tsvg -o vcf_body.svg vcf_body.dot`)
=#

const vcf_metainfo_actions = Dict(
    :metainfo_key      => :(record.key = (mark1:p-1) - offset),
    :metainfo_val      => :(record.val = (mark2:p-1) - offset; record.dict = data[mark2] == UInt8('<')),
    :metainfo_dict_key => :(push!(record.dictkey, (mark1:p-1) - offset)),
    :metainfo_dict_val => :(push!(record.dictval, (mark1:p-1) - offset)),
    :metainfo          => :(),
    :header_sampleID   => :(push!(reader.header.sampleID, String(data[mark1:p-1]))),
    :vcfheader         => :(found_header = true; @escape),
    :countline         => :(linenum += 1),
    :anchor            => :(),
    :mark              => :(mark1 = p),
    :mark2             => :(mark2 = p))

eval(Bio.ReaderHelper.generate_index_function(VCFMetaInfo, vcf_metainfo_machine, vcf_metainfo_actions))

const vcf_record_actions = Dict(
    :record_chrom        => :(record.chrom = (mark:p-1) - offset),
    :record_pos          => :(record.pos = (mark:p-1) - offset),
    :record_id           => :(push!(record.id, (mark:p-1) - offset)),
    :record_ref          => :(record.ref = (mark:p-1) - offset),
    :record_alt          => :(push!(record.alt, (mark:p-1) - offset)),
    :record_qual         => :(record.qual = (mark:p-1) - offset),
    :record_filter       => :(push!(record.filter, (mark:p-1) - offset)),
    :record_info_key     => :(push!(record.infokey, (mark:p-1) - offset)),
    :record_format       => :(push!(record.format, (mark:p-1) - offset)),
    :record_genotype     => :(push!(record.genotype, UnitRange{Int}[])),
    :record_genotype_elm => :(push!(record.genotype[end], (mark:p-1) - offset)),
    :record              => :(),
    :countline           => :(linenum += 1),
    :anchor              => :(),
    :mark                => :(mark = p))

eval(Bio.ReaderHelper.generate_index_function(VCFRecord, vcf_record_machine, vcf_record_actions))

const vcf_header_actions = merge(vcf_metainfo_actions, Dict(
    :metainfo => quote
        record.data = data[Bio.ReaderHelper.upanchor!(stream):p-1]
        record.filled = true
        push!(reader.header.metainfo, record)
        record = VCFMetaInfo()
    end,
    :anchor    => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))

function readheader!(reader::VCFReader)
    _readheader!(reader, reader.state)
end

@eval function _readheader!(reader::VCFReader, state::Bio.Ragel.State)
    stream = state.stream
    Bio.ReaderHelper.ensure_margin!(stream)
    cs = state.cs
    linenum = state.linenum
    data = stream.buffer
    p = stream.position
    p_end = stream.available
    p_eof = -1
    offset = mark1 = mark2 = 0
    found_header = false
    record = VCFMetaInfo()

    while true
        $(Automa.generate_exec_code(vcf_header_machine, actions=vcf_header_actions, code=:table))

        @assert cs != 0
        state.cs = cs
        state.finished = cs == 0
        state.linenum = linenum
        stream.position = p

        if cs < 0
            error("VCF file format error on line ", linenum)
        elseif found_header
            Bio.ReaderHelper.upanchor!(stream)
            break
        #elseif cs == 0
        #    throw(EOFError())
        elseif p > p_eof ≥ 0
            error("incomplete VCF input on line ", linenum)
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

const vcf_body_actions = merge(vcf_record_actions, Dict(
    :record => :(found_record = true; @escape),
    :anchor => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))

eval(Bio.ReaderHelper.generate_read_function(VCFReader, vcf_body_machine, vcf_body_actions))
