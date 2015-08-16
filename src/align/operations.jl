# =====================
# Alignment operations
# =====================


# Single operations are encoded in one byte each
# ----------------------------------------------
bitstype 8 Operation


# Conversion to and from integers
# -------------------------------

convert(::Type{Operation}, num::Uint8) = box(Operation, unbox(Uint8, num))
convert(::Type{Uint8}, op::Operation) = box(Uint8, unbox(Operation, op))

convert{T<:Unsigned}(::Type{Operation}, unint::T) = convert(Operation, convert(Uint8, unint))
convert{T<:Unsigned}(::Type{T}, op::Operation) = convert(T, convert(Uint8, op))


# Operation encoding definitions
# ------------------------------

#= Meanings of the Operations
-----------------------------

+--------------+----------------+--------------------------------------------+
| Operation    | CIGAR Rep      | Desription                                 |
+--------------+----------------+--------------------------------------------+
| OP_GAP       | Not Applicale. | Denotes a plain stretch of gap characters. |
| OP_MM        | CIGAR 'M'      | Match/Mismatch.                            |
| OP_N         | CIGAR 'N'      |                                            |
| OP_MATCH     | CIGAR '='      | Match.                                     |
| OP_MISMATCH  | CIGAR 'X'      | Mismatch.                                  |
| OP_SCLIP     | CIGAR 'S'      | Soft clipping.                             |
| OP_HCLIP     | CIGAR 'H'      | Hard clipping.                             |
| OP_INSERT    | CIGAR 'I'      | Insertion.                                 |
| OP_DELETE    | CIGAR 'D'      | Deletion.                                  |
| OP_PAD       | CIGAR 'P'      | Padding.                                   |
| OP_INVALID   | Invalid        | Invalid operation.                         |
+--------------+----------------+--------------------------------------------+

=#

for op = [
  (:OP_GAP, 0), (:OP_MM, 1), (:OP_N, 2), (:OP_MATCH, 3),
  (:OP_MISMATCH, 4), (:OP_SCLIP, 5), (:OP_HCLIP, 6), (:OP_INSERT, 7),
  (:OP_DELETE, 8), (:OP_PAD, 9), (:OP_INVALID, 255)
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
  OP_INVALID, OP_MATCH, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_DELETE, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_HCLIP, OP_INSERT, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_MM, OP_N, OP_INVALID,
  OP_PAD, OP_INVALID, OP_INVALID, OP_SCLIP, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_MISMATCH, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_DELETE, OP_INVALID, OP_INVALID, OP_INVALID, OP_HCLIP,
  OP_INSERT, OP_INVALID, OP_INVALID, OP_INVALID, OP_MM,
  OP_N, OP_INVALID, OP_PAD, OP_INVALID, OP_INVALID,
  OP_SCLIP, OP_INVALID, OP_INVALID, OP_INVALID, OP_INVALID,
  OP_MISMATCH
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
