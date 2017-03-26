# VCF Reader
# ==========

type Reader <: Bio.IO.AbstractReader
    state::Bio.Ragel.State
    header::Header

    function Reader(input::BufferedStreams.BufferedInputStream)
        reader = new(Bio.Ragel.State(vcf_header_machine.start_state, input), Header())
        readheader!(reader)
        reader.state.cs = vcf_body_machine.start_state
        return reader
    end
end

"""
    VCF.Reader(input::IO)

Create a data reader of the VCF file format.

# Arguments
* `input`: data source
"""
function Reader(input::IO)
    return Reader(BufferedStreams.BufferedInputStream(input))
end

function Base.eltype(::Type{Reader})
    return Record
end

function Bio.IO.stream(reader::Reader)
    return reader.state.stream
end

"""
    header(reader::VCF.Reader)::VCF.Header

Get the header of `reader`.
"""
function header(reader::Reader)
    return reader.header
end

function Bio.header(reader::Reader)
    return header(reader)
end

# VCF v4.3
Base.info("compiling VCF")
const vcf_metainfo_machine, vcf_record_machine, vcf_header_machine, vcf_body_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    delim(x, sep) = opt(cat(x, rep(cat(sep, x))))

    # The 'fileformat' field is required and must be the first line.
    fileformat = let
        key = cat("fileformat")
        key.actions[:enter] = [:mark1]
        key.actions[:exit]  = [:metainfo_tag]

        version = re"[!-~]+"
        version.actions[:enter] = [:mark2]
        version.actions[:exit]  = [:metainfo_val]

        cat("##", key, '=', version)
    end
    fileformat.actions[:enter] = [:anchor]
    fileformat.actions[:exit]  = [:metainfo]

    # All kinds of meta-information line after 'fileformat' are handled here.
    metainfo = let
        tag = re"[0-9A-Za-z_]+"
        tag.actions[:enter] = [:mark1]
        tag.actions[:exit]  = [:metainfo_tag]

        str = re"[ -;=-~][ -~]*"  # does not starts with '<'
        str.actions[:enter] = [:mark2]
        str.actions[:exit]  = [:metainfo_val]

        dict = let
            dictkey = re"[0-9A-Za-z_]+"
            dictkey.actions[:enter] = [:mark1]
            dictkey.actions[:exit]  = [:metainfo_dict_key]

            dictval = let
                quoted   = cat('"', rep(alt(re"[ !#-[\]-~]", "\\\"", "\\\\")), '"')
                unquoted = rep(re"[ -~]" \ re"[\",>]")
                alt(quoted, unquoted)
            end
            dictval.actions[:enter] = [:mark1]
            dictval.actions[:exit]  = [:metainfo_dict_val]

            cat('<', delim(cat(dictkey, '=', dictval), ','), '>')
        end
        dict.actions[:enter] = [:mark2]
        dict.actions[:exit]  = [:metainfo_val]

        cat("##", tag, '=', alt(str, dict))
    end
    metainfo.actions[:enter] = [:anchor]
    metainfo.actions[:exit]  = [:metainfo]

    # The header line.
    header = let
        sampleID = re"[ -~]+"
        sampleID.actions[:enter] = [:mark1]
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

const vcf_metainfo_actions = Dict(
    :metainfo_tag      => :(record.tag = (mark1:p-1) - offset),
    :metainfo_val      => :(record.val = (mark2:p-1) - offset; record.dict = data[mark2] == UInt8('<')),
    :metainfo_dict_key => :(push!(record.dictkey, (mark1:p-1) - offset)),
    :metainfo_dict_val => :(push!(record.dictval, (mark1:p-1) - offset)),
    :metainfo          => quote
        Bio.ReaderHelper.resize_and_copy!(record.data, data, offset+1:p-1)
        record.filled = (offset+1:p-1) - offset
    end,
    :anchor            => :(),
    :mark1             => :(mark1 = p),
    :mark2             => :(mark2 = p))
eval(
    Bio.ReaderHelper.generate_index_function(
        MetaInfo,
        vcf_metainfo_machine,
        :(mark1 = mark2 = offset = 0),
        vcf_metainfo_actions))
eval(
    Bio.ReaderHelper.generate_readheader_function(
        Reader,
        MetaInfo,
        vcf_header_machine,
        :(mark1 = mark2 = offset = 0),
        merge(vcf_metainfo_actions, Dict(
            :metainfo => quote
                Bio.ReaderHelper.resize_and_copy!(record.data, data, Bio.ReaderHelper.upanchor!(stream):p-1)
                record.filled = (offset+1:p-1) - offset
                @assert isfilled(record)
                push!(reader.header.metainfo, record)
                Bio.ReaderHelper.ensure_margin!(stream)
                record = MetaInfo()
            end,
            :header_sampleID => :(push!(reader.header.sampleID, String(data[mark1:p-1]))),
            :vcfheader => :(finish_header = true; @escape),
            :countline => :(linenum += 1),
            :anchor    => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))))

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
    :record              => quote
        Bio.ReaderHelper.resize_and_copy!(record.data, data, 1:p-1)
        record.filled = (offset+1:p-1) - offset
    end,
    :anchor              => :(),
    :mark                => :(mark = p))
eval(
    Bio.ReaderHelper.generate_index_function(
        Record,
        vcf_record_machine,
        :(mark = offset = 0),
        vcf_record_actions))
eval(
    Bio.ReaderHelper.generate_read_function(
        Reader,
        vcf_body_machine,
        :(mark = offset = 0),
        merge(vcf_record_actions, Dict(
            :record    => quote
                Bio.ReaderHelper.resize_and_copy!(record.data, data, Bio.ReaderHelper.upanchor!(stream):p-1)
                record.filled = (offset+1:p-1) - offset
                found_record = true
                @escape
            end,
            :countline => :(linenum += 1),
            :anchor    => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))))
