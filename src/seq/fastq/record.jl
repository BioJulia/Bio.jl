# FASTQ Record
# ============

type Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # indexes
    identifier::UnitRange{Int}
    description::UnitRange{Int}
    sequence::UnitRange{Int}
    quality::UnitRange{Int}
end

"""
    FASTQ.Record()

Create an unfilled FASTQ record.
"""
function Record()
    return Record(UInt8[], 1:0, 1:0, 1:0, 1:0, 1:0)
end

"""
    FASTQ.Record(data::Vector{UInt8})

Create a FASTQ record object from `data`.

This function verifies and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(data, 1:0, 1:0, 1:0, 1:0, 1:0)
    index!(record)
    return record
end

"""
    FASTQ.Record(str::AbstractString)

Create a FASTQ record object from `str`.

This function verifies and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, convert(Vector{UInt8}, str))
end

"""
    FASTQ.Record(identifier, sequence, quality; offset=33)

Create a FASTQ record from `identifier`, `sequence` and `quality`.
"""
function Record(identifier::AbstractString, sequence, quality::Vector; offset=33)
    return Record(identifier, nothing, sequence, quality, offset=offset)
end

"""
    FASTQ.Record(identifier, description, sequence, quality; offset=33)

Create a FASTQ record from `identifier`, `description`, `sequence` and `quality`.
"""
function Record(identifier::AbstractString, description::Union{AbstractString,Void}, sequence, quality::Vector; offset=33)
    if length(sequence) != length(quality)
        throw(ArgumentError("the length of sequence doesn't match the length of quality"))
    end
    buf = IOBuffer()
    print(buf, '@', identifier)
    if description != nothing
        print(buf, ' ', description)
    end
    print(buf, '\n')
    print(buf, sequence, '\n')
    print(buf, "+\n")
    ascii_quality = convert(Vector{UInt8}, quality + offset)
    write(buf, ascii_quality, '\n')
    return Record(takebuf_array(buf))
end

function Base.convert(::Type{Record}, str::AbstractString)
    return convert(Record, convert(Vector{UInt8}, str))
end

function Base.copy(record::Record)
    return Record(
        record.data[record.filled],
        record.filled,
        record.identifier,
        record.description,
        record.sequence,
        record.quality)
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   identifier: ", hasidentifier(record) ? identifier(record) : "<missing>")
        println(io, "  description: ", hasdescription(record) ? description(record) : "<missing>")
        println(io, "     sequence: ", hassequence(record) ? sequence(String, record) : "<missing>")
          print(io, "      quality: ", hasquality(record) ? quality(record) : "<missing>")
    else
        print(io, " <not filled>")
    end
end

function initialize!(record::Record)
    record.filled = 1:0
    record.identifier = 1:0
    record.description = 1:0
    record.sequence = 1:0
    record.quality = 1:0
    return record
end

function Bio.isfilled(record::Record)
    return !isempty(record.filled)
end


# Accessor functions
# ------------------

"""
    identifier(record::Record)::String

Get the sequence identifier of `record`.
"""
function identifier(record::Record)::String
    checkfilled(record)
    if !hasidentifier(record)
        missingerror(:identifier)
    end
    return String(record.data[record.identifier])
end

function hasidentifier(record::Record)
    return isfilled(record)
end

"""
    description(record::Record)::String

Get the description of `record`.
"""
function description(record::Record)::String
    checkfilled(record)
    if !hasdescription(record)
        missingerror(:description)
    end
    return String(record.data[record.description])
end

function hasdescription(record)
    return isfilled(record) && record.description != 1:0
end

"""
    sequence(::Type{S}, record::Record, [part::UnitRange{Int}])::S

Get the sequence of `record`.

`S` can be either a subtype of `Bio.Seq.Sequence` or `String`.
If `part` argument is given, it returns the specified part of the sequence.
"""
function sequence{S<:Bio.Seq.Sequence}(::Type{S}, record::Record, part::UnitRange{Int}=1:endof(record.sequence))::S
    checkfilled(record)
    if !hassequence(record)
        missingerror(:sequence)
    end
    seqpart = record.sequence[part]
    return S(record.data, first(seqpart), last(seqpart))
end

function sequence(::Type{String}, record::Record, part::UnitRange{Int}=1:endof(record.sequence))::String
    checkfilled(record)
    if !hassequence(record)
        missingerror(:sequence)
    end
    return String(record.data[record.sequence[part]])
end

"""
    sequence(record::Record, [part::UnitRange{Int}])::Bio.Seq.DNASequence

Get the sequence of `record`.
"""
function sequence(record::Record, part::UnitRange{Int}=1:endof(record.sequence))::Bio.Seq.DNASequence
    checkfilled(record)
    return sequence(Bio.Seq.DNASequence, record, part)
end

function hassequence(record::Record)
    # zero-length sequence may exist
    return isfilled(record)
end

"""
    quality(record::Record, [offset::Integer=33, [part::UnitRange]])::Vector{UInt8}

Get the base quality of `record`.
"""
function quality(record::Record, offset::Integer=33, part::UnitRange{Int}=1:endof(record.quality))::Vector{UInt8}
    checkfilled(record)
    quality = record.data[record.quality[part]]
    for i in 1:endof(part)
        # TODO: Checked arithmetic?
        @inbounds quality[i] -= offset
    end
    return quality
end

"""
    quality(record::Record, encoding_name::Symbol, [part::UnitRange])::Vector{UInt8}

Get the base quality of `record` by decoding with `encoding_name`.

The `encoding_name` can be either `:sanger`, `:solexa`, `:illumina13`, `:illumina15`, or `:illumina18`.
"""
function quality(record::Record, encoding_name::Symbol, part::UnitRange{Int}=1:endof(record.quality))::Vector{UInt8}
    checkfilled(record)
    encoding = (
        encoding_name == :sanger     ?     SANGER_QUAL_ENCODING :
        encoding_name == :solexa     ?     SOLEXA_QUAL_ENCODING :
        encoding_name == :illumina13 ? ILLUMINA13_QUAL_ENCODING :
        encoding_name == :illumina15 ? ILLUMINA15_QUAL_ENCODING :
        encoding_name == :illumina18 ? ILLUMINA18_QUAL_ENCODING :
        throw(ArgumentError("quality encoding ':$(encoding_name)' is not supported")))
    quality = Vector{UInt8}(length(part))
    if !isempty(part)
        qpart = record.quality[part]
        check_quality_string(encoding, record.data, first(qpart), last(qpart))
        decode_quality_string!(encoding, record.data, quality, first(qpart), last(qpart))
    end
    return quality
end

function hasquality(record::Record)
    return isfilled(record)
end

function Bio.seqname(record::Record)
    return identifier(record)
end

function Bio.hasseqname(record::Record)
    return hasidentifier(record)
end

function Bio.sequence(record::Record)
    return sequence(record)
end

function Bio.hassequence(record::Record)
    return hassequence(record)
end


# Helper functions
# ----------------

function checkfilled(record)
    if !isfilled(record)
        throw(ArgumentError("unfilled FASTQ record"))
    end
end
