# =====================
# Alignment operations
# =====================


# Single operations are encoded in one byte each
# ----------------------------------------------
bitstype 8 Operation


# Conversion to and from integers
# -------------------------------

convert(::Type{Operation}, num::UInt8) = box(Operation, unbox(UInt8, num))
convert(::Type{UInt8}, op::Operation) = box(UInt8, unbox(Operation, op))

convert{T<:Unsigned}(::Type{Operation}, unint::T) = convert(Operation, convert(UInt8, unint))
convert{T<:Unsigned}(::Type{T}, op::Operation) = convert(T, convert(UInt8, op))


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


function convert(::Type{Operation}, c::Char)
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

function convert(::Type{Char}, op::Operation)
    @assert op != OP_INVALID error("Alignment operation is not valid.")
    return op_to_char[convert(UInt8, op) + 1]
end

function show(io::IO, op::Operation)
    write(io, convert(Char, op))
end
