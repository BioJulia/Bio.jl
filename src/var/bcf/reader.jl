# BCF Reader
# ==========

immutable Reader{T<:IO} <: Bio.IO.AbstractReader
    version::Tuple{UInt8,UInt8}  # (major, minor)
    header::Bio.Var.VCF.Header
    stream::BGZFStreams.BGZFStream{T}
end

"""
    BCF.Reader(input::IO)

Create a data reader of the BCF file format.

# Arguments
* `input`: data source
"""
function Reader(input::IO)
    stream = BGZFStreams.BGZFStream(input)

    # magic bytes and BCF version
    B = read(stream, UInt8)
    C = read(stream, UInt8)
    F = read(stream, UInt8)
    major = read(stream, UInt8)
    minor = read(stream, UInt8)
    if (B, C, F) != map(UInt8, ('B', 'C', 'F'))
        error("not a BCF file")
    elseif (major, minor) != (0x02, 0x02)
        error("unsupported BCF version")
    end

    # NOTE: This isn't specified in the BCF specs, but seems to be a practice at least in htslib.
    # See: https://github.com/samtools/hts-specs/issues/138
    l_header = read(stream, Int32)
    data = read(stream, l_header)

    # parse VCF header
    vcfreader = Bio.Var.VCF.Reader(BufferedStreams.BufferedInputStream(data))

    return Reader((major, minor), vcfreader.header, stream)
end

function Base.eltype{T}(::Type{Reader{T}})
    return Record
end

function Bio.IO.stream(reader::Reader)
    return reader.stream
end

"""
    header(reader::BCF.Reader)::VCF.Header

Get the header of `reader`.
"""
function header(reader::Reader)
    return reader.header
end

function Bio.header(reader::Reader)
    return header(reader)
end

function Base.read!(reader::Reader, record::Record)
    sharedlen = read(reader.stream, UInt32)
    indivlen = read(reader.stream, UInt32)
    datalen = sharedlen + indivlen
    resize!(record.data, datalen)
    unsafe_read(reader.stream, pointer(record.data), datalen)
    record.filled = 1:datalen
    record.sharedlen = sharedlen
    record.indivlen = indivlen
    return record
end
