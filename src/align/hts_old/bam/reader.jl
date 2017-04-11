# BAM Reader
# ==========

"""
    BAMReader(input::IO; index=nothing)

Create a data reader of the BAM file format.

# Arguments
* `input`: data source
* `index=nothing`: filepath to a random access index (currently *bai* is Supported)
"""
type BAMReader{T} <: Bio.IO.AbstractReader
    stream::BGZFStreams.BGZFStream{T}
    header::SAMHeader
    start_offset::BGZFStreams.VirtualOffset
    refseqnames::Vector{String}
    refseqlens::Vector{Int}
    index::Nullable{BAI}
end

function BAMReader(input::IO; index=nothing)
    if isa(index, AbstractString)
        index = BAI(index)
    else
        if index != nothing
            error("unrecognizable index argument")
        end
    end
    reader = init_bam_reader(input)
    reader.index = index
    return reader
end

function Base.eltype{T}(::Type{BAMReader{T}})
    return BAMRecord
end

function Bio.IO.stream(reader::BAMReader)
    return reader.stream
end

function Base.show(io::IO, reader::BAMReader)
    println(io, summary(reader), ":")
    println(io, "  header keys: ", join(keys(reader.header), ", "))
      print(io, "  number of contigs: ", length(reader.refseqnames))
end

function header(reader::BAMReader, fillSQ::Bool=false)
    if fillSQ
        if haskey(reader.header, "SQ")
            # TODO: check consistency
        else
            records = Dict{String,Any}[]
            for (name, len) in zip(reader.refseqnames, reader.refseqlens)
                push!(records, Dict("SN" => name, "LN" => len))
            end
            reader.header["SQ"] = records
        end
    end
    return reader.header
end

function Base.seek(reader::BAMReader, voffset::BGZFStreams.VirtualOffset)
    seek(reader.stream, voffset)
end

function Base.seekstart(reader::BAMReader)
    seek(reader.stream, reader.start_offset)
end

# Initialize a BAM reader by reading the header section.
function init_bam_reader(input::IO)
    stream = BGZFStreams.BGZFStream(input)

    # magic bytes
    B = read(stream, UInt8)
    A = read(stream, UInt8)
    M = read(stream, UInt8)
    x = read(stream, UInt8)
    if B != UInt8('B') || A != UInt8('A') || M != UInt8('M') || x != 0x01
        error("input was not a valid BAM file")
    end

    # SAM header
    textlen = read(stream, Int32)
    samheader = parse_samheader(read(stream, UInt8, textlen))

    # reference sequences
    refseqnames = String[]
    refseqlens = Int[]
    n_refs = read(stream, Int32)
    for _ in 1:n_refs
        namelen = read(stream, Int32)
        data = read(stream, UInt8, namelen)
        seqname = unsafe_string(pointer(data))
        seqlen = read(stream, Int32)
        push!(refseqnames, seqname)
        push!(refseqlens, seqlen)
    end

    reader = BAMReader(
        stream,
        samheader,
        isa(input, Pipe) ? BGZFStreams.VirtualOffset(0, 0) : BGZFStreams.virtualoffset(stream),
        refseqnames,
        refseqlens,
        Nullable{BAI}())

    return reader
end

function Base.read!(reader::BAMReader, record::BAMRecord)
    unsafe_read(
        reader.stream,
        pointer_from_objref(record),
        BAM_FIXED_FIELDS_BYTES)
    dsize = data_size(record)
    if length(record.data) < dsize
        resize!(record.data, dsize)
    end
    unsafe_read(reader.stream, pointer(record.data), dsize)
    record.refseqnames = reader.refseqnames
    return record
end

function Base.start(reader::BAMReader)
    return BAMRecord()
end

function Base.done(reader::BAMReader, rec)
    return eof(reader)
end

function Base.next(reader::BAMReader, rec)
    read!(reader, rec)
    return copy(rec), rec
end
