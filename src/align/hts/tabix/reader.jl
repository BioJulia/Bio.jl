# Tabix Reader
# ============

function Base.read(filename::AbstractString, ::Type{Tabix})
    return open(input -> read(input, Tabix), filename)
end

function Base.read(input_::IO, ::Type{Tabix})
    input = BGZFStream(input_)

    # magic
    T = read(input, UInt8)
    B = read(input, UInt8)
    I = read(input, UInt8)
    x = read(input, UInt8)
    if T != UInt8('T') || B != UInt8('B') || I != UInt8('I') || x != 0x01
        error("invalid tabix magic bytes")
    end

    n_refs = read(input, Int32)
    format = read(input, Int32)
    col_seq = read(input, Int32)
    col_beg = read(input, Int32)
    col_end = read(input, Int32)
    meta = read(input, Int32)
    skip = read(input, Int32)

    l_nm = read(input, Int32)
    data = read(input, UInt8, l_nm)
    names = split(String(data), '\0', keep=false)

    indexes = read_indexes(input, n_refs)

    if !eof(input)
        n_no_coor = Nullable(read(input, UInt64))
    else
        n_no_coor = Nullable{UInt64}()
    end

    return Tabix(
        format,
        (col_seq, col_beg, col_end),
        meta,
        skip,
        names,
        indexes,
        n_no_coor)
end
