# The zlib Compression
# ====================

function compress!(dst::Vector{UInt8}, src::Vector{UInt8})::UInt64
    sz = Ref(Culong(length(dst)))
    code = ccall(
        (:compress, Libz.zlib),
        Cint,
        (Ptr{Void}, Ref{Culong}, Ptr{Void}, Culong),
        dst, sz, src, sizeof(src))
    if code != Libz.Z_OK
        Libz.zerror(code)
    end
    return sz[]
end

function uncompress!(dst::Vector{UInt8}, src::Vector{UInt8})::UInt64
    sz = Ref(Culong(length(dst)))
    code = ccall(
        (:uncompress, Libz.zlib),
        Cint,
        (Ptr{Void}, Ref{Culong}, Ptr{Void}, Culong),
        dst, sz, src, sizeof(src))
    if code != Libz.Z_OK
        Libz.zerror(code)
    end
    return sz[]
end
