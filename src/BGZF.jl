
module BGZF

export BGZFSource

using BufferedStreams
using Libz

#import Base: eof, readbytes!


# compressed and decompressed blocks are <= this in BGZF
#const BGZF_MAX_BLOCK_SIZE = Int(0x10000)  # FIXME: 64 * 1024
const BGZF_MAX_BLOCK_SIZE = 64 * 1024


"""
General error thrown when data did not conform the BGZF standard.
"""
immutable MalformedBGZFData <: Exception
end


immutable BGZFSource{T <: IO}
    input::T
    zstream::Base.RefValue{Libz.ZStream}

    # FIXME: remove compressed_block_ptr
    # space to read the next compressed block
    compressed_block::Vector{UInt8}
    compressed_block_ptr::Ptr{UInt8}

    # FIXME: remove decompressed_block_ptr
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
    zstream = Libz.init_inflate_zstream(true)
    # FIXME: Vector{UInt8}(...)
    compressed_block = Vector{UInt8}(BGZF_MAX_BLOCK_SIZE)
    decompressed_block = Vector{UInt8}(BGZF_MAX_BLOCK_SIZE)
    return BGZFSource(input, zstream,
                      compressed_block, pointer(compressed_block),
                      decompressed_block, pointer(decompressed_block),
                      Ref(0), Ref(0), Ref(false))
end


# FIXME: Base.eof
function Base.eof(source::BGZFSource)
    return source.eof[]
end


#=
"""
Read the next compressed BGZF block into `output`.
"""
function read_bgzf_block!(source::BGZFSource)
    input = source.input

    # read header up to xlen
    id1 = read(input, UInt8)
    id2 = read(input, UInt8)
    cm  = read(input, UInt8)
    flg = read(input, UInt8)

    # 0x1f = 31, 0x8b = 139
    if id1 != 0x1f || id2 != 0x8b || cm != 0x08 || flg != 0x04
        throw(MalformedBGZFData)
    end

    # FIXME: explicitly read MTIME, XFL, and OS
    seekforward(input, 6)
    xlen = Int(read(input, UInt16))
    if xlen < 6
        throw(MalformedBGZFData)
    end

    # read extra subfields
    nb = seekforward(input, 4)
    bsize = Int(read(input, UInt16))
    nb += 2
    nb += seekforward(input, xlen - 6)
    if nb != xlen
        throw(MalformedBGZFData)
    end

    # read the rest of the bgzf block
    output = source.compressed_block
    remaining_block_size = bsize - xlen - 11
    nb = readbytes!(input, output, remaining_block_size)
    if nb != remaining_block_size
        throw(MalformedBGZFData)
    end

    # size of uncompressed (input) data
    isize = (Int(output[nb]) << 24) | (Int(output[nb - 1]) << 16) |
            (Int(output[nb - 2]) << 8) | (Int(output[nb - 3]))
    return remaining_block_size - 8, isize
end
=#


"""
Decompress the next BGZF block into source.decompressed_block.
"""
function decompress_block(source::BGZFSource)
    #bsize, isize = read_bgzf_block!(source)
    bsize = read_bgzf_block!(source)

    zstream = getindex(source.zstream)
    #zstream.next_out = source.decompressed_block_ptr
    #zstream.avail_out = isize
    #zstream.next_in = source.compressed_block_ptr
    #zstream.avail_in = bsize
    zstream.next_out = pointer(source.decompressed_block)
    zstream.avail_out = BGZF_MAX_BLOCK_SIZE
    zstream.next_in = pointer(source.compressed_block)
    zstream.avail_in = bsize

    ret = ccall((:inflate, Libz._zlib), Cint, (Ptr{Libz.ZStream}, Cint),
                source.zstream, Libz.Z_FINISH)

    if ret != Libz.Z_STREAM_END || zstream.avail_in != 0 # || zstream.avail_out != 0
        error("Failed to decompress a BGZF block (zlib error $(ret))")
    end
   n_avail = BGZF_MAX_BLOCK_SIZE - zstream.avail_out

    ret = ccall((:inflateReset, Libz._zlib), Cint, (Ptr{Libz.ZStream},), source.zstream)
    if ret != Libz.Z_OK
        error("Unable to reset zlib stream.")
    end

    source.bytes_consumed[] = 0
    #source.bytes_available[] = isize
    source.bytes_available[] = n_avail
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
            source.eof[] = true
        else
            throw(MalformedBGZFData)
        end
    end
    return bsize
end

function BufferedStreams.readbytes!(source::BGZFSource, buffer::Vector{UInt8}, from::Int, to::Int)
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

#function BufferedStreams.readbytes!(source::BGZFSource, buffer::Vector{UInt8}, from::Int, to::Int)
#    
#end

if VERSION < v"0.5-"
    function readbytesto!(s::IO, buffer::Vector{UInt8},
                            from::Integer, nb::Integer)
        ptr = pointer(buffer, from)
        buffrom = pointer_to_array(ptr, length(buffer) - from + 1, false)
        return readbytes!(s, buffrom, nb)
    end
else
    function readbytesto!(s::IO, buffer::Vector{UInt8},
                            from::Integer, nb::Integer)
        return unsafe_read(io, pointer(buffer, from), nb)
    end
end

end # module BGZF
