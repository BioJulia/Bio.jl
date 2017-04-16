# FASTQ Reader
# ============

immutable Reader <: Bio.IO.AbstractReader
    state::Bio.Ragel.State
    seq_transform::Nullable{Function}

    function Reader(input::BufferedInputStream, seq_transform)
        return new(Bio.Ragel.State(file_machine.start_state, input), seq_transform)
    end
end

"""
    FASTQ.Reader(input::IO; fill_ambiguous=nothing)

Create a data reader of the FASTQ file format.

# Arguments
* `input`: data source
* `fill_ambiguous=nothing`: fill ambiguous symbols with the given symbol
"""
function Reader(input::IO; fill_ambiguous=nothing)
    if fill_ambiguous === nothing
        seq_transform = nothing
    else
        seq_transform = generate_fill_ambiguous(fill_ambiguous)
    end
    return Reader(BufferedInputStream(input), seq_transform)
end

function Base.eltype(::Type{Reader})
    return Record
end

function Bio.IO.stream(reader::Reader)
    return reader.state.stream
end

function generate_fill_ambiguous(symbol::BioSymbols.DNA)
    certain = map(UInt8, ('A', 'C', 'G', 'T', 'a', 'c', 'g', 't'))
    # return transform function
    return function (data, range)
        fill = convert(UInt8, convert(Char, symbol))
        for i in range
            if data[i] ∉ certain
                data[i] = fill
            end
        end
        return data
    end
end

# NOTE: This does not support line-wraps within sequence and quality.
info("compiling FASTQ")
const record_machine, file_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    rep1 = Automa.RegExp.rep1
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    any = Automa.RegExp.any
    space = Automa.RegExp.space

    hspace = re"[ \t\v]"

    header1 = let
        identifier = rep(any() \ space())
        identifier.actions[:enter] = [:mark]
        identifier.actions[:exit]  = [:header1_identifier]

        description = cat(any() \ hspace, re"[^\r\n]*")
        description.actions[:enter] = [:mark]
        description.actions[:exit]  = [:header1_description]

        cat('@', identifier, opt(cat(rep1(hspace), description)))
    end

    sequence = re"[A-z]*"
    sequence.actions[:enter] = [:mark]
    sequence.actions[:exit]  = [:sequence]

    header2 = let
        identifier = rep1(any() \ space())
        identifier.actions[:enter] = [:mark]
        identifier.actions[:exit]  = [:header2_identifier]

        description = cat(any() \ hspace, re"[^\r\n]*")
        description.actions[:enter] = [:mark]
        description.actions[:exit]  = [:header2_description]

        cat('+', opt(cat(identifier, opt(cat(rep1(hspace), description)))))
    end

    quality = re"[!-~]*"
    quality.actions[:enter] = [:mark]
    quality.actions[:exit]  = [:quality]

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(opt('\r'), lf)
    end

    record′ = cat(header1, newline, sequence, newline, header2, newline, quality)
    record′.actions[:enter] = [:anchor]
    record′.actions[:exit]  = [:record]
    record = cat(record′, newline)

    file = rep(record)

    return map(Automa.compile, (record, file))
end)()

#=
write("fastq.dot", Automa.machine2dot(file_machine))
run(`dot -Tsvg -o fastq.svg fastq.dot`)
=#

function check_identical(data1, range1, data2, range2)
    if length(range1) != length(range2) ||
       memcmp(pointer(data1, first(range1)), pointer(data2, first(range2)), length(range1)) != 0
       error("sequence and quality have non-matching header")
    end
end

function memcmp(p1::Ptr, p2::Ptr, len::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, len) % Int
end

const record_actions = Dict(
    :header1_identifier  => :(record.identifier  = (mark:p-1)),
    :header1_description => :(record.description = (mark:p-1)),
    :header2_identifier  => :(check_identical(record.data, mark:p-1, record.data, record.identifier)),
    :header2_description => :(check_identical(record.data, mark:p-1, record.data, record.description)),
    :sequence => :(record.sequence = (mark:p-1)),
    :quality  => :(record.quality  = (mark:p-1)),
    :record   => :(record.filled   = 1:p-1),
    :anchor => :(),
    :mark   => :(mark = p),
    :countline => :())
eval(
    Bio.ReaderHelper.generate_index_function(
        Record,
        record_machine,
        :(mark = 0),
        record_actions))
eval(
    Bio.ReaderHelper.generate_read_function(
        Reader,
        file_machine,
        :(mark = offset = 0),
        merge(record_actions, Dict(
            :header1_identifier  => :(record.identifier  = (mark:p-1) - stream.anchor + 1),
            :header1_description => :(record.description = (mark:p-1) - stream.anchor + 1),
            :header2_identifier  => :(check_identical(data, mark:p-1, data, (record.identifier) + stream.anchor - 1)),
            :header2_description => :(check_identical(data, mark:p-1, data, (record.description) + stream.anchor - 1)),
            :sequence            => :(record.sequence    = (mark:p-1) - stream.anchor + 1),
            :quality             => :(record.quality     = (mark:p-1) - stream.anchor + 1),
            :record => quote
                if length(record.sequence) != length(record.quality)
                    error("the length of sequence does not match the length of quality")
                end
                Bio.ReaderHelper.resize_and_copy!(record.data, data, Bio.ReaderHelper.upanchor!(stream):p-1)
                record.filled = (offset+1:p-1) - offset
                if !isnull(reader.seq_transform)
                    get(reader.seq_transform)(record.data, record.sequence)
                end
                found_record = true
                @escape
            end,
            :countline => :(linenum += 1),
            :anchor => :(Bio.ReaderHelper.anchor!(stream, p); offset = p - 1)))))
