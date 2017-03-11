# FASTA Record
# ============

type Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # indexes
    identifier::UnitRange{Int}
    description::UnitRange{Int}
    sequence::UnitRange{Int}
end

"""
    FASTA.Record()

Create an unfilled FASTA record.
"""
function Record()
    return Record(UInt8[], 1:0, 1:0, 1:0, 1:0)
end

"""
    FASTA.Record(data::Vector{UInt8})

Create a FASTA record object from `data`.
This function verifies and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(data, 1:0, 1:0, 1:0, 1:0)
    index!(record)
    return record
end

"""
    FASTA.Record(str::AbstractString)

Create a FASTA record object from `str`.
This function verifies and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return Record(convert(Vector{UInt8}, str))
end

function Base.convert(::Type{String}, record::Record)
    return String(record.data[datarange(record)])
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = datarange(record1)
        r2 = datarange(record2)
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function Base.copy(record::Record)
    return Record(record.data[record.filled], record.identifier, record.description, record.sequence)
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled), length(record.filled)))
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
          print(io, "     sequence: ", hassequence(record) ? sequence(String, record) : "<missing>")
    else
        print(io, " <not filled>")
    end
end

function initialize!(record::Record)
    record.filled = 1:0
    record.identifier = 1:0
    record.description = 1:0
    record.sequence = 1:0
    return record
end

function isfilled(record::Record)
    return !isempty(record.filled)
end

function datarange(record::Record)
    return record.filled
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
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

function hasidentifier(record)
    return isfilled(record) && !isempty(record.identifier)
end

function Bio.seqname(record::Record)
    return identifier(record)
end

function Bio.hasseqname(record::Record)
    return hasidentifier(record)
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
    return isfilled(record) && !isempty(record.description)
end

"""
    sequence(::Type{S}, record::Record)::S

Get the sequence of `record`.
`S` can be either a subtype of `Bio.Seq.Sequence` or `String`.
"""
function sequence{S<:Bio.Seq.Sequence}(::Type{S}, record::Record)::S
    checkfilled(record)
    if !hassequence(record)
        missingerror(:sequence)
    end
    return convert(S, record.data[record.sequence])
end

function sequence(::Type{String}, record::Record)::String
    checkfilled(record)
    if !hassequence(record)
        missingerror(:sequence)
    end
    return String(record.data[record.sequence])
end

function hassequence(record::Record)
    # zero-length sequence may exist
    return isfilled(record)
end

function checkfilled(record)
    if !isfilled(record)
        throw(ArgumentError("unfilled FASTA record"))
    end
end
