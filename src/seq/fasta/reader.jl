# FASTA Reader
# ============

type Reader <: Bio.IO.AbstractReader
    state::Bio.Ragel.State
    index::Nullable{Index}

    function Reader(input::BufferedStreams.BufferedInputStream, index)
        return new(Bio.Ragel.State(file_machine.start_state, input), index)
    end
end

"""
    FASTAReader(input::IO; index=nothing)

Create a data reader of the FASTA file format.

# Arguments
* `input`: data source
* `index=nothing`: filepath to a random access index (currently *fai* is supported)
"""
function Reader(input::IO; index=nothing)
    return Reader(BufferedStreams.BufferedInputStream(input), index)
end

function Base.eltype(::Type{Reader})
    return Record
end

function Bio.IO.stream(reader::Reader)
    return reader.state.stream
end

info("compiling FASTA")
const record_machine, file_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    rep1 = Automa.RegExp.rep1
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    any = Automa.RegExp.any
    space = Automa.RegExp.space

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(opt('\r'), lf)
    end

    hspace = re"[ \t\v]"
    whitespace = space() | newline

    identifier = rep1(any() \ space())
    identifier.actions[:enter] = [:mark]
    identifier.actions[:exit]  = [:identifier]

    description = cat(any() \ hspace, re"[^\r\n]*")
    description.actions[:enter] = [:mark]
    description.actions[:exit]  = [:description]

    header = cat('>', identifier, opt(cat(rep1(hspace), description)))
    header.actions[:enter] = [:anchor]
    header.actions[:exit]  = [:header]

    letters = rep1(any() \ space() \ cat('>'))
    letters.actions[:enter] = [:anchor, :mark]
    letters.actions[:exit]  = [:letters]

    sequence = opt(cat(letters, rep(cat(rep1(newline), letters))))
    sequence.actions[:enter] = [:sequence_start]

    record = cat(header, rep1(newline), sequence, rep1(newline))
    record.actions[:enter] = [:mark1]
    record.actions[:exit]  = [:record]

    file = cat(rep(newline), rep(record))

    return map(Automa.compile, (record, file))
end)()

#=
write("fasta.dot", Automa.machine2dot(file_machine))
run(`dot -Tsvg -o fasta.svg fasta.dot`)
=#

const record_actions = Dict(
    :identifier  => :(record.identifier  = (mark:p-1)),
    :description => :(record.description = (mark:p-1)),
    :header      => quote
        @assert record.data === data
        copy!(record.data, 1, record.data, 1, p - mark1)
        copied += p - mark1
        if p ≤ endof(data)
            data[p] = UInt8('\n')
            copied += 1
        end
    end,
    :letters => quote
        len = p - mark
        copy!(record.data, copied + 1, record.data, mark, p - mark)
        copied += len
        record.sequence = first(record.sequence):last(record.sequence)+len
    end,
    :sequence_start => :(record.sequence = copied+1:copied),
    :record => :(record.filled = 1:copied),
    :anchor => :(),
    :mark => :(mark = p),
    :mark1 => :(mark1 = p),
    :countline => :(#= linenum += 1 =#))
eval(
    Bio.ReaderHelper.generate_index_function(
        Record,
        record_machine,
        record_actions))
eval(
    Bio.ReaderHelper.generate_read_function(
        Reader,
        file_machine,
        merge(record_actions, Dict(
            :identifier  => :(record.identifier  = (mark:p-1) - stream.anchor + 1),
            :description => :(record.description = (mark:p-1) - stream.anchor + 1),
            :header => quote
                range = Bio.ReaderHelper.upanchor!(stream):p-1
                Bio.ReaderHelper.resize_and_copy!(record.data, copied + 1, data, range)
                copied += length(range)
                if copied + 1 ≤ endof(record.data)
                    record.data[copied] = UInt8('\n')
                else
                    push!(record.data, UInt8('\n'))
                end
                copied += 1
            end,
            :sequence_start => :(record.sequence = copied+1:copied),
            :letters => quote
                range = Bio.ReaderHelper.upanchor!(stream):p-1
                Bio.ReaderHelper.resize_and_copy!(record.data, copied + 1, data, range)
                record.sequence = first(record.sequence):last(record.sequence)+length(range)
                copied += length(range)
            end,
            :record => quote
                record.filled = 1:copied
                found_record = true
                @escape
            end,
            :countline => :(linenum += 1),
            :anchor => :(Bio.ReaderHelper.anchor!(stream, p))))))
