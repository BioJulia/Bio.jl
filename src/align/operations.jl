# Alignment Operations
# ====================
#
# Alignment operation type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# An alignment operation type
bitstype 8 Operation


# Conversion to and from integers
# -------------------------------

Base.convert(::Type{Operation}, op::UInt8) = reinterpret(Operation, op)
Base.convert(::Type{UInt8}, op::Operation) = reinterpret(UInt8, op)

Base.convert{T<:Unsigned}(::Type{Operation}, op::T) = Operation(UInt8(op))
Base.convert{T<:Unsigned}(::Type{T}, op::Operation) = T(UInt8(op))


# Operation encoding definitions
# ------------------------------

const OP_MATCH        = convert(Operation, 0x00) # M
const OP_INSERT       = convert(Operation, 0x01) # I
const OP_DELETE       = convert(Operation, 0x02) # D
const OP_SKIP         = convert(Operation, 0x03) # N
const OP_SOFT_CLIP    = convert(Operation, 0x04) # S
const OP_HARD_CLIP    = convert(Operation, 0x05) # H
const OP_PAD          = convert(Operation, 0x06) # P
const OP_SEQ_MATCH    = convert(Operation, 0x07) # =
const OP_SEQ_MISMATCH = convert(Operation, 0x08) # X
const OP_BACK         = convert(Operation, 0x09) # B
const OP_START        = convert(Operation, 0x0a) # 0 (non-standard)
const OP_INVALID      = convert(Operation, 0xff)

const OP_MAX_VALID = OP_START

# classify operations
function ismatchop(op::Operation)
    return op == OP_MATCH || op == OP_SEQ_MATCH || op == OP_SEQ_MISMATCH
end

function isinsertop(op::Operation)
    return op == OP_INSERT || op == OP_SOFT_CLIP || op == OP_HARD_CLIP
end

function isdeleteop(op::Operation)
    return op == OP_DELETE || op == OP_SKIP
end


# Conversion from characters to operations
# ----------------------------------------

# Lookup table for conversion from Char to Operation
const char_to_op = [
    OP_SEQ_MATCH, OP_INVALID, OP_INVALID,   OP_INVALID,
    OP_INVALID,   OP_BACK,    OP_INVALID,   OP_DELETE,
    OP_INVALID,   OP_INVALID, OP_INVALID,   OP_HARD_CLIP,
    OP_INSERT,    OP_INVALID, OP_INVALID,   OP_INVALID,
    OP_MATCH,     OP_SKIP,    OP_INVALID,   OP_PAD,
    OP_INVALID,   OP_INVALID, OP_SOFT_CLIP, OP_INVALID,
    OP_INVALID,   OP_INVALID, OP_INVALID,   OP_SEQ_MISMATCH ]

function Base.convert(::Type{Operation}, c::Char)
    @inbounds op = '=' <= c <= 'X' ? char_to_op[c - '=' + 1] :
                          c == '0' ? OP_START                : OP_INVALID
    if op == OP_INVALID
        error("Invalid alignment operation '$c'")
    end
    return op
end


# Conversion from characters to operations
# ----------------------------------------

const op_to_char = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '0']

function Base.convert(::Type{Char}, op::Operation)
    @assert op != OP_INVALID error("Alignment operation is not valid.")
    return op_to_char[convert(UInt8, op) + 1]
end

function Base.show(io::IO, op::Operation)
    write(io, convert(Char, op))
    return
end
