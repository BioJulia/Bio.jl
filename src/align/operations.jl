# =====================
# Alignment operations
# =====================

# Operation types are encoded in one byte each
# ---------------------------------------------
bitstype 8 OpKind


# Conversion to and from integers
# -------------------------------

convert(::Type{OpKind}, num::Uint8) = box(OpKind, unbox(Uint8, num))
convert(::Type{Uint8}, op::OpKind) = box(Uint8, unbox(OpKind, op))

convert{T<:Unsigned}(::Type{OpKind}, unint::T) = convert(OpKind, convert(Uint8, unint))
convert{T<:Unsigned}(::Type{T}, op::OpKind) = convert(T, convert(Uint8, op))


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

    @eval const $(op[1]) = convert(OpKind, Uint8($(op[2])))

end


# Conversion from characters to operation type
# ---------------------------------------------

# Lookup table for conversion from Char to Operation
const char_to_op_kind = [
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

function convert(::Type{OpKind}, c::Char)
    @inbounds op = '-' <= c <= 'x' ? char_to_op_kind[c - '-' + 1] : OP_INVALID
    if op == OP_INVALID
        error("$(c) is not a valid alignment operation.")
    end
    return op
end


# Conversion from characters to operation type
# ---------------------------------------------

const op_kind_to_char = ['-', 'M', 'N', '=', 'X', 'S', 'H', 'I', 'D', 'P']

function convert(::Type{Char}, op::OpKind)
    @assert op != OP_INVALID error("Alignment operation is not valid.")
    @inbounds ch = op_kind_to_char[convert(Uint8, op) + 1]
    return ch
end

function show(io::IO, op::OpKind)
    write(io, convert(Char, op))
end


immutable Operation
    Sort::OpKind
    Size::Int
end

function Operation(kind::Char, size::Int)
    return Operation(OpKind(kind), size)
end

copy(x::Operation) = Operation(x.Sort, x.Size)

function convert(::Type{Operation}, str::String)
    return Operation(str[end], parse(Int, str[1:end-1]))
end

function convert(::Type{String}, operation::Operation)
    return string(operation.Size, Char(operation.Sort))
end

function show(io::IO, operation::Operation)
    write(io, convert(String, operation))
end
