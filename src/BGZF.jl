
module BGZF

export BGZFSource

using BufferedStreams
using Libz

import Base: eof, readbytes!


# compressed and decompressed blocks are <= this in BGZF
const BGZF_MAX_BLOCK_SIZE = Int(0x10000)


"""
General error thrown when data did not conform the BGZF standard.
"""
immutable MalformedBGZFData <: Exception
end


immutable BGZFSource{T <: IO}
    input::T
    zstream::Base.RefValue{Libz.ZStream}

    # space to read the next compressed block
    compressed_block::Vector{UInt8}
    compressed_block_ptr::Ptr{UInt8}

    # space to decompress the block
    decompressed_block::Vector{UInt8}
    decompressed_block_ptr::Ptr{UInt8}

    # number of bytes available in decompressed_block
    bytes_available::Base.RefValue{Int}

    # number of bytes from decompressed_block that have been consumed
    bytes_consumed::Base.RefValue{Int}

    # true if there are no more bytes to consume
    eof::Base.RefValue{Bool}
end


function BGZFSource(input::IO)
    zstream = Libz.init_inflate_zstream(false)
    compressed_block = Array(UInt8, BGZF_MAX_BLOCK_SIZE)
    decompressed_block = Array(UInt8, BGZF_MAX_BLOCK_SIZE)
    return BGZFSource(input, zstream,
                      compressed_block, pointer(compressed_block),
                      decompressed_block, pointer(decompressed_block),
                      Ref(0), Ref(0), Ref(false))
end


function eof(source::BGZFSource)
    return source.eof[]
end


"""
Read the next compressed BGZF block into `output`.
"""
function read_bgzf_block!(source::BGZFSource)
    input = source.input
    output = source.compressed_block

    # read header up to xlen
    p = 1
    id1 = output[p] = read(input, UInt8); p += 1
    id2 = output[p] = read(input, UInt8); p += 1
    cm  = output[p] = read(input, UInt8); p += 1
    flg = output[p] = read(input, UInt8); p += 1

    if id1 != 0x1f || id2 != 0x8b || cm != 0x08 || flg != 0x04
        throw(MalformedBGZFData)
    end

    for _ in 1:8
        output[p] = read(input, UInt8); p += 1
    end

    xlen = (Int(output[p - 1]) << 8) | Int(output[p - 2])
    if xlen < 6
        throw(MalformedBGZFData)
    end

    # read extra subfields
    nb = readbytes!(input, output, p, p + xlen - 1)
    if nb != xlen
        throw(MalformedBGZFData)
    end
    bsize = (Int(output[p + 5]) << 8 | Int(output[p + 4]))
    p += xlen

    # read the rest of the bgzf block
    remaining_block_size = bsize - xlen - 11
    nb = readbytes!(input, output, p, p + remaining_block_size - 1)
    p += nb
    if nb != remaining_block_size
        throw(MalformedBGZFData)
    end

    # size of uncompressed (input) data
    isize = (Int(output[p - 1]) << 24) | (Int(output[p - 2]) << 16) |
            (Int(output[p - 3]) << 8) | (Int(output[p - 4]))
    return bsize + 1, isize
end


"""
Decompress the next BGZF block into source.decompressed_block.
"""
function decompress_block(source::BGZFSource)
    bsize, isize = read_bgzf_block!(source)

    zstream = getindex(source.zstream)
    zstream.next_out = source.decompressed_block_ptr
    zstream.avail_out = isize
    zstream.next_in = source.compressed_block_ptr
    zstream.avail_in = bsize

    ret = ccall((:inflate, Libz._zlib), Cint, (Ptr{Libz.ZStream}, Cint),
                source.zstream, Libz.Z_FINISH)

   if ret != Libz.Z_STREAM_END || zstream.avail_in != 0 || zstream.avail_out != 0
       error("Failed to decompress a BGZF block (zlib error $(ret))")
   end

    ret = ccall((:inflateReset, Libz._zlib), Cint, (Ptr{Libz.ZStream},), source.zstream)
    if ret != Libz.Z_OK
        error("Unable to reset zlib stream.")
    end

    source.bytes_consumed[] = 0
    source.bytes_available[] = isize
end


function readbytes!(source::BGZFSource, buffer::Vector{UInt8}, from::Int, to::Int)
    from0 = from
    while to - from + 1 > 0
        available = source.bytes_available[] - source.bytes_consumed[]
        nb = min(to - from + 1, available)
        Base.unsafe_copy!(pointer(buffer, from),
                          pointer(source.decompressed_block,
                                  1 + source.bytes_consumed[]), nb)
        from += nb
        source.bytes_consumed[] += nb

        if source.bytes_consumed[] == source.bytes_available[]
            if eof(source.input)
                source.eof[] = true
                break
            else
                decompress_block(source)
            end
        end
    end

    return from - from0
end


end # module BGZF

