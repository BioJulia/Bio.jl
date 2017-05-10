# Memory-level Operations
# =======================

module Mem

function copy(dst, src, n::Integer)
    return copy(pointer(dst), pointer(src), n)
end

function copy(dst::Ptr, src::Ptr, n::Integer)
    return ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Csize_t), dst, src, n)
end

end
