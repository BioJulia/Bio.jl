# Assign flags for sequence types so we can maintain a set of compatible
# alphabets in an integer.

"""
Type representing an alphabet

An `Alphabet` value holds a set of alphabets compatible with a sequence.
Usually this is just one.
"""
bitstype 16 Alphabet

function convert(::Type{Alphabet}, nt::UInt16)
    return box(Alphabet, unbox(UInt16, nt))
end


function convert(::Type{UInt16}, nt::Alphabet)
    return box(UInt16, unbox(Alphabet, nt))
end


function (|)(a::Alphabet, b::Alphabet)
    return convert(Alphabet, convert(UInt16, a) | convert(UInt16, b))
end


function (&)(a::Alphabet, b::Alphabet)
    return convert(Alphabet, convert(UInt16, a) & convert(UInt16, b))
end

# for safe module precompilation
Base.hash(a::Alphabet) = hash(convert(UInt16, a))

"`Alphabet` value indicating no compatible alphabets."
const EMPTY_ALPHABET = convert(Alphabet, UInt16(0))

"DNA alphabet"
const DNA_ALPHABET   = convert(Alphabet, UInt16(0b0001))

"RNA alphabet"
const RNA_ALPHABET   = convert(Alphabet, UInt16(0b0010))

"amino acid alphabet"
const AA_ALPHABET    = convert(Alphabet, UInt16(0b0100))

"`Alphabet` value indicating that all known alphabets are compatible"
const ALL_ALPHABETS =
    DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET


const alphabet_type = Dict{Alphabet, Type}(
    DNA_ALPHABET => DNASequence,
    RNA_ALPHABET => RNASequence,
    AA_ALPHABET  => AminoAcidSequence
)


# When a sequence has multiple compatible alphabets, we choose the first
# compatible alphabet in this list.
const preferred_sequence_alphabets = [
    DNA_ALPHABET, RNA_ALPHABET, AA_ALPHABET
]


# Lookup table mapping a character in 'A':'z' to an integer representing the set
# of alphabets that character is compatible with.
const compatible_alphabets = [
    # A                                          B
      DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # C                                          D
      DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # E            F            G
      AA_ALPHABET, AA_ALPHABET, DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET,
    # H            I            J            K            L
      AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # M            N                                          O
      AA_ALPHABET, DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # P            Q            R            S
      AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # T                           U                           V
      DNA_ALPHABET | AA_ALPHABET, RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # W,           X            Y            Z
      AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET,

    # [               \               ]               ^
      EMPTY_ALPHABET, EMPTY_ALPHABET, EMPTY_ALPHABET, EMPTY_ALPHABET,
    #  _              `
      EMPTY_ALPHABET, EMPTY_ALPHABET,

    # a                                          b
      DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # c                                          d
      DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # e            f            g
      AA_ALPHABET, AA_ALPHABET, DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET,
    # h            i            j            k            l
      AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # m            n                                          o
      AA_ALPHABET, DNA_ALPHABET | RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # p            q            r            s
      AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET,
    # t                           u                           v
      DNA_ALPHABET | AA_ALPHABET, RNA_ALPHABET | AA_ALPHABET, AA_ALPHABET,
    # w            x            y            z
      AA_ALPHABET, AA_ALPHABET, AA_ALPHABET, AA_ALPHABET
]


"""
Infer the sequence type by inspecting a string.

### Arguments
   * `data`: sequence data in a string
   * `start`: first position to consider in data
   * `stop`: last position to consider in data
   * `default`: if there are multiple compatible alphabets, default
             to this one if it's compatible.

### Returns
A type T to which the string data can be converted.
"""
function infer_alphabet(data::Vector{UInt8}, start, stop, default)
    alphabets = ALL_ALPHABETS
    if start > stop
        return default
    end

    if start < 1 ||  stop > length(data)
        throw(BoundsError())
    end

    @inbounds for i in start:stop
        c = data[i]
        if UInt8('A') <= c <= UInt8('z')
            alphabets &= compatible_alphabets[c - UInt8('A') + 1]
        else
            error("Character $(c) is not compatible with any sequence type.")
        end
    end

    if count_ones(convert(UInt16, alphabets)) == 0
        error("String is not compatible with any known sequence type.")
    elseif alphabets & default != EMPTY_ALPHABET
        return default
    else
        for alphabet in preferred_sequence_alphabets
            if alphabet & alphabets != EMPTY_ALPHABET
                return alphabet
            end
        end
    end
    return default
end
