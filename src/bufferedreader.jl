import Base: read, read!, seek, position, eof, nb_available, eof


const BUFFERED_READER_INITIAL_BUF_SIZE = 1000000

type BufferedReader <: IO
    # Input source. If null, it's assumed that buffer contains
    # all the input.
    input::Nullable{IO}
    input_position::UInt

    # Should we close input when finished
    input_owned::Bool

    # Buffered data.
    buffer::Vector{UInt8}

    # Position of the last position in buffer that's filled.
    buffer_end::UInt

    # An index into buffer that will be updated and stay valid upon refill.
    mark::UInt
end

function BufferedReader(data::Vector{UInt8}, memory_map::Bool=false, len=length(data))
    if memory_map
        error("Parser must be given a file name in order to memory map.")
    end
    return BufferedReader(Nullable{IO}(), 0, false, data, len, 0)
end

function BufferedReader(input::IO, memory_map::Bool=false)
    if memory_map
        error("Parser must be given a file name in order to memory map.")
    end
    return BufferedReader(Nullable{IO}(input), 1, false,
                          Array(UInt8, BUFFERED_READER_INITIAL_BUF_SIZE), 0, 0)
end

function BufferedReader(filename::String, memory_map::Bool=false)
    if memory_map
        data = Mmap.mmap(open(filename), Vector{UInt8}, (filesize(filename),))
        return BufferedReader(Nullable{IO}(), 0, false, data, length(data), 0)
    else
        return BufferedReader(Nullable{IO}(open(filename)), 1, true,
                              Array(UInt8, BUFFERED_READER_INITIAL_BUF_SIZE),
                              0, 0)
    end
end


"""
Refill a buffer, keeping some portion of it.

The currently set mark will be shifted to the beginning of the buffer, the the
rest of the buffer will be filled by reading bytes from the input.

This is useful if a parser is the middle of matching something when the buffer
runs out

# Args
  `reader`: A BufferedReader.

# Modifies
  reader.buffer is refilled, and indexes are updated to the correct positions
  in the refilled buffer

# Returns
A (newpos, num_bytes_read), where newpos is the index within the newly
refilled buffer that corresponds to the end of the buffer prior to refilling. If
the buffer was entirely refilled and nothing was kept this will be 0.
`num_bytes_read` is the number of bytes read from the input source.
"""
function fillbuffer!(reader::BufferedReader)
    if isnull(reader.input)
        return UInt(0)
    end

    buflen = length(reader.buffer)
    keeplen = UInt(0)
    if reader.mark > 0
        keeplen = reader.buffer_end - reader.mark + 1
        if keeplen == buflen
            buflen = 2 * buflen
            resize!(reader.buffer, buflen)
        end
        copy!(reader.buffer, 1, reader.buffer, reader.mark, keeplen)
        reader.mark = 1
    end

    nb = readchunk!(get(reader.input), reader.buffer, keeplen + 1, buflen)
    reader.buffer_end = keeplen + nb

    if nb == 0 && reader.input_owned
        close(get(reader.input))
    end

    reader.input_position += nb

    return nb
end


function readchunk!(source::IO, dest::Vector{UInt8}, dest_start::Integer,
                    dest_stop::Integer)
    i = dest_start
    while i <= dest_stop && !eof(source)
        @inbounds dest[i] = read(source, UInt8)
        i += 1
    end
    return UInt(i - dest_start)
end

function readchunk!(source::IOStream, dest::Vector{UInt8}, dest_start::Integer,
                    dest_stop::Integer)
    len = dest_stop - dest_start  + 1
    return ccall(:ios_readall, UInt, (Ptr{Void}, Ptr{Void}, UInt), source.ios,
                 pointer(dest, dest_start), len)
end

function mark!(reader::BufferedReader, mark::Integer)
    reader.mark = mark
end

function unmark!(reader::BufferedReader)
    mark = reader.mark
    reader.mark = 0
    return mark
end

# IO Interface
# ------------

# We implement an IO interface on top of BufferedReader by treading buffer.mark
# as the current position within the stream.

"""
Seek to the given position in the input stream. Unlike Base.seek, this is 1-based.
"""
function seek(reader::BufferedReader, pos::Integer)
    if isnull(reader.input)
        if pos <= length(reader.buffer)
            reader.mark = pos
        else
            throw(BoundsError)
        end
    else
        if applicable(seek, get(reader.input), pos)
            if reader.input_position - reader.buffer_end <= pos < reader.input_position
                reader.mark = pos - (reader.input_position - reader.buffer_end) + 1
            else
                seek(get(reader.input), pos - 1)
                reader.mark = 1
                reader.buffer_end = 0
            end
        else
            error("Can't seek in input stream.")
            # TODO: Allow seeking forwards
        end
    end
end

function eof(reader::BufferedReader)
    return reader.mark > reader.buffer_end &&
        (isnull(reader.input) || eof(get(reader.input)))
end

"""
Current position in the stream. 1-based, unlike Base.seek.
"""
function position(reader::BufferedReader)
    @assert reader.mark != 0
    # TODO: better double check this shit
    return reader.input_position - (reader.buffer_end - reader.mark)
end

function read(reader::BufferedReader, ::Type{UInt8})
    if reader.mark == 0
        reader.mark = 1
    end

    if reader.mark > reader.buffer_end
        reader.mark = 0
        nb = fillbuffer!(reader)
        if nb == 0
            error("Unexpected end of input")
        end
        reader.mark = 1
    end

    b = reader.buffer[reader.mark]
    reader.mark += 1
    return b
end

function nb_available(reader::BufferedReader)
    if reader.buffer_end >= reader.mark
        return reader.buffer_end - reader.mark + 1
    else
        return 0
    end
end

function eof(reader::BufferedReader)
    if reader.mark > reader.buffer_end
        return fillbuffer!(reader) == 0
    else
        return false
    end
end
