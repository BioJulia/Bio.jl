# BGZF
# ====
#
# The BGZF compression format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module BGZF

export BGZFSource, BGZFSink

using BufferedStreams
using Libz

# compressed and decompressed blocks are <= this in BGZF
const BGZF_MAX_BLOCK_SIZE = 64 * 1024


"""
General error thrown when data did not conform the BGZF standard.
"""
immutable MalformedBGZFData <: Exception
end


type BGZFSource{T<:IO}
    input::T
    zstream::Base.RefValue{Libz.ZStream}

    # space to read the next compressed block
    compressed_block::Vector{UInt8}

    # space to decompress the block
    decompressed_block::Vector{UInt8}

    # number of bytes available in decompressed_block
    bytes_available::Int

    # number of bytes from decompressed_block that have been consumed
    bytes_consumed::Int

    # true if there are no more bytes to consume
    eof::Bool
end

function BGZFSource(input::IO)
    zstream = Libz.init_inflate_zstream(true)
    compressed_block = Vector{UInt8}(BGZF_MAX_BLOCK_SIZE)
    decompressed_block = Vector{UInt8}(BGZF_MAX_BLOCK_SIZE)
    return BGZFSource(input, zstream,
                      compressed_block,
                      decompressed_block,
                      0, 0, false)
end

Base.eof(source::BGZFSource) = source.eof

# decompress the next BGZF block into `source.decompressed_block`.
function decompress_block(source::BGZFSource)
    bsize = read_bgzf_block!(source)

    zstream = getindex(source.zstream)
    zstream.next_out = pointer(source.decompressed_block)
    zstream.avail_out = BGZF_MAX_BLOCK_SIZE
    zstream.next_in = pointer(source.compressed_block)
    zstream.avail_in = bsize

    ret = ccall((:inflate, Libz._zlib), Cint, (Ptr{Libz.ZStream}, Cint),
                source.zstream, Libz.Z_FINISH)

    if ret != Libz.Z_STREAM_END || zstream.avail_in != 0 # || zstream.avail_out != 0
        error("failed to decompress a BGZF block (zlib error $(ret))")
    end
   n_avail = BGZF_MAX_BLOCK_SIZE - zstream.avail_out

    ret = ccall((:inflateReset, Libz._zlib), Cint, (Ptr{Libz.ZStream},), source.zstream)
    if ret != Libz.Z_OK
        error("failed to reset zlib stream.")
    end

    source.bytes_consumed = 0
    source.bytes_available = n_avail
    return
end

function read_bgzf_block!(source::BGZFSource)
    input = source.input
    block = source.compressed_block

    # +---+---+---+---+---+---+---+---+---+---+
    # |ID1|ID2|CM |FLG|     MTIME     |XFL|OS | (more-->)
    # +---+---+---+---+---+---+---+---+---+---+
    n = readbytes!(input, block, 10)
    if n != 10 || block[1] != 0x1f || block[2] != 0x8b || block[3] != 0x08 || block[4] != 0x04
        throw(MalformedBGZFData)
    end

    # +---+---+=================================+
    # | XLEN  |...XLEN bytes of "extra field"...| (more-->)
    # +---+---+=================================+
    n = readbytesto!(input, block, 11, 2)
    if n != 2
        throw(MalformedBGZFData)
    end
    xlen = UInt16(block[11]) | UInt16(block[12]) << 8
    n = readbytesto!(input, block, 13, xlen)
    if n != xlen
        throw(MalformedBGZFData)
    end
    bsize::UInt16 = 0
    pos = 12
    while pos < 12 + xlen
        si1 = block[pos+1]
        si2 = block[pos+2]
        slen = UInt16(block[pos+3]) | UInt16(block[pos+4]) << 8
        if si1 == 0x42 || si2 == 0x43
            if slen != 2
                throw(MalformedBGZFData)
            end
            bsize = (UInt16(block[pos+5]) | UInt16(block[pos+6]) << 8) + 1
        end
        # skip this field
        pos += 4 + slen
    end
    if bsize == 0
        # extra subfied of the BGZF file format is not found
        throw(MalformedBGZFData)
    end

    # +=======================+---+---+---+---+---+---+---+---+
    # |...compressed blocks...|     CRC32     |     ISIZE     |
    # +=======================+---+---+---+---+---+---+---+---+
    readbytesto!(input, block, 13 + xlen, bsize - (12 + xlen))

    if eof(input)
        # the last block must be an end-of-file marker
        if bsize == 28
            source.eof = true
        else
            throw(MalformedBGZFData)
        end
    end
    return bsize
end

function BufferedStreams.readbytes!(source::BGZFSource, buffer::Vector{UInt8}, from::Int, to::Int)
    from0 = from
    while to - from + 1 > 0
        available = source.bytes_available - source.bytes_consumed
        nb = min(to - from + 1, available)
        Base.unsafe_copy!(pointer(buffer, from),
                          pointer(source.decompressed_block,
                                  1 + source.bytes_consumed), nb)
        from += nb
        source.bytes_consumed += nb

        if source.bytes_consumed == source.bytes_available
            if eof(source.input)
                source.eof = true
                break
            else
                decompress_block(source)
            end
        end
    end

    return from - from0
end

function readbytesto!(s::IO, buffer::Vector{UInt8},
                        from::Integer, nb::Integer)
    ptr = pointer(buffer, from)
    buffrom = pointer_to_array(ptr, length(buffer) - from + 1, false)
    return readbytes!(s, buffrom, nb)
end


type BGZFSink{T<:IO}
    output::T
    zstream::Base.RefValue{Libz.ZStream}

    # space to read the next compressed block
    compressed_block::Vector{UInt8}

    # space to decompress the block
    decompressed_block::Vector{UInt8}
end

function BGZFSink(output::IO)
    zstream = Libz.init_deflate_stream(
        true, Libz.Z_DEFAULT_COMPRESSION, 8, Libz.Z_DEFAULT_STRATEGY)
    compressed_block = Vector{UInt8}(BGZF_MAX_BLOCK_SIZE)
    decompressed_block = Vector{UInt8}(BGZF_MAX_BLOCK_SIZE)
    return BGZFSink(output, zstream, compressed_block, decompressed_block)
end

function BufferedStreams.writebytes(sink::BGZFSink, buffer::Vector{UInt8}, n::Int, eof::Bool)
    # ...
end

end  # module BGZF
