# =====================
# Alignment operations
# =====================


# Single operations are encoded in one byte each
# ----------------------------------------------
bitstype 8 Operation


# Conversion to and from integers
# -------------------------------

convert(::Type{Operation}, num::Uint8)  = box(Operation, unbox(Uint8, num))
convert(::Type{Uint8}, op::Operation) = box(Uint8, unbox(Operation, op))

convert{T<:Unsigned}(::Type{Operation}, unint::T) = convert(Operation, convert(Uint8, unint))
convert{T<:Unsigned}(::Type{T}, op::Operation) = box(T, Base.zext_int(T, unbox(Operation, op)))


# Operation encoding definitions
# ------------------------------

const OP_GAP = convert(Operation, Uint8(0))         # Denotes a plain gap region
const OP_CIGAR_M = convert(Operation, Uint8(1))     # CIGAR M | Match/Mismatch
const OP_CIGAR_N = convert(Operation, Uint8(2))     # CIGAR N
const OP_CIGAR_EQ = convert(Operation, Uint8(3))    # CIGAR = | Match
const OP_CIGAR_X = convert(Operation, Uint8(4))     # CIGAR X | Mismatch
const OP_CIGAR_S = convert(Operation, Uint8(5))     # CIGAR S | Soft clipping
const OP_CIGAR_H = convert(Operation, Uint8(6))     # CIGAR H | Hard clipping
const OP_CIGAR_I = convert(Operation, Uint8(7))     # CIGAR I | Insertion
const OP_CIGAR_D = convert(Operation, Uint8(8))     # CIGAR D | Deletion
const OP_CIGAR_P = convert(Operation, Uint8(9))     # CIGAR P | Padding
const OP_INVALID = convert(Operation, Uint8(255))   # Invalid Operation


# Conversion from characters to operations
# ----------------------------------------

# Lookup table for conversion from Char to Operation
const char_to_op = [
  OP_GAP, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID, 
  OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_INVALID, OP_CIGAR_EQ, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_CIGAR_D, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_CIGAR_H, OP_CIGAR_I, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_CIGAR_M, OP_CIGAR_N, OP_INVALID,
  OP_CIGAR_P, OP_INVALID, OP_INVALID, OP_CIGAR_S, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_CIGAR_X, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_CIGAR_D, OP_INVALID, OP_INVALID, OP_INVALID, OP_CIGAR_H,
  OP_CIGAR_I, OP_INVALID, OP_INVALID, OP_INVALID, OP_CIGAR_M,
  OP_CIGAR_N, OP_INVALID, OP_CIGAR_P, OP_INVALID, OP_INVALID,
  OP_CIGAR_S, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_CIGAR_X 
]


function convert(::Type{Operation}, c::Char)
    @inbounds op = '-' <= c <= 'x' ? char_to_op[c - '-' + 1] : OP_INVALID
    @assert op != OP_INVALID error("$(c) is not a valid alignment operation.")
    return op
end


# Conversion from characters to operations
# ----------------------------------------

const op_to_char = ['-', 'M', 'N', '=', 'X', 'S', 'H', 'I', 'D', 'P']

function convert(::Type{Char}, op::Operation)
  @assert op != OP_INVALID error("Alignment operation is not valid.")
  @inbounds ch = op_to_char[convert(Uint8, op) + 1]
  return ch
end

function show(io::IO, op::Operation)
    write(io, convert(Char, op))
end
