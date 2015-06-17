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
convert{T<:Unsigned}(::Type{T}, op::Operation) = convert(T, convert(Uint8, op))


# Operation encoding definitions
# ------------------------------

#= Meanings of the Operations
-----------------------------

OP_GAP      | '-'       | Denotes a plain gap region
OP_CIGAR_M  | CIGAR 'M' | Match/Mismatch
OP_CIGAR_N  | CIGAR 'N' |
OP_CIGAR_EQ | CIGAR '=' | Match
OP_CIGAR_X  | CIGAR 'X' | Mismatch
OP_CIGAR_S  | CIGAR 'S' | Soft clipping
OP_CIGAR_H  | CIGAR 'H' | Hard clipping
OP_CIGAR_I  | CIGAR 'I' | Insertion
OP_CIGAR_D  | CIGAR 'D' | Deletion
OP_CIGAR_P  | CIGAR 'P' | Padding
OP_INVALID  | Invalid   | Invalid operation.

=#

for op = [
  (:OP_GAP, 0), (:OP_CIGAR_M, 1), (:OP_CIGAR_N, 2), (:OP_CIGAR_EQ, 3),
  (:OP_CIGAR_X, 4), (:OP_CIGAR_S, 5), (:OP_CIGAR_H, 6), (:OP_CIGAR_I, 7),
  (:OP_CIGAR_D, 8), (:OP_CIGAR_P, 9), (:OP_INVALID, 255)
  ]

  @eval const $(op[1]) = convert(Operation, Uint8($(op[2])))

end


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
