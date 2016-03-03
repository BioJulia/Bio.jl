# Alphabet
# ========
#
# Subtypes of `Alphabet` represent a domain of biological characters. For
# example, `DNAAlphabet{2}` has a domain of unambiguous nucleotides (i.e. A, C,
# G, and T). These types are used for parameterizing biological sequences and so
# on. A pair of encoder and decoder is associated with an alphabet, which maps
# values between binary and Julia-level representation.

"""
Alphabet of biological characters.
"""
abstract Alphabet

"""
The number of bits to represent the alphabet.
"""
function bitsof end

# Creating an alphabet requires quite a few type and method definitions for each alphabet.
# This macro makes it simpler to define an alphabet with less boilerplate code.
macro alphabet(name, bits, element_type, charset)
    code = quote end

    # Name should be a symbol or expression.
    type_name = typeof(name) == Symbol ? name : Symbol(name)
    # Bits should be a single integer or array of integers.
    multi_encodings = false
    if typeof(bits) <: Integer
        bit_encoding = bits
    else
        if bits.head == :vect
            bit_encoding = bits.args
            multi_encodings = true
        end
    end

    # Eltype should be a single value, or an array.
    element_type =
        typeof(element_type) == Symbol ? Symbol[element_type] : convert(Array{Symbol, 1}, element_type.args)
    multi_eltypes =
        typeof(element_type) == Array && length(element_type) > 1 ? true : false

    # Alphabet should be a single range/function/array,
    # or an array of functions/ranges/arrays.
    if typeof(charset) != Expr
        error("Invalid options provided as charset")
    end
    multi_charsets = false
    charset_head = charset.head
    if charset_head == :vect
        charset = charset.args
        types = DataType[typeof(c) for c in charset]
        if length(unique(types)) == 1 && all(types .!= Expr)
            multi_charsets = false
            charset = Array[charset]
        elseif length(charset) == length(bit_encoding)
            multi_charsets = true
        else
            error("Invalid options provided as charset")
        end
    elseif charset_head == :call || charset_head == :(:)
        charset = Expr[charset]
    end

    # Generate the code for the alphabet.
    typetodef = multi_encodings ? :($(type_name){n}) : :($type_name)
    push!(code.args, esc(:(immutable $(typetodef) <: Alphabet end)))
    just_type = :(::Type{$type_name})
    for i in bit_encoding
        t = multi_encodings ? :(::Type{$type_name{$i}}) : just_type
        push!(code.args, esc(:(function bitsof($t)
                               $i
                               end)))
    end
    for i in 1:length(element_type)
        sn = multi_eltypes ? bit_encoding[i] : :(n)
        fun = sn == :(n) && multi_encodings ? :(Base.eltype{n}) : :(Base.eltype)
        argument = multi_encodings ? :(::Type{$type_name{$sn}}) : just_type
        push!(code.args, esc(:(function $fun($argument)
                               $(element_type[i])
                               end)))
    end
    for i in 1:length(charset)
        sn = multi_charsets ? bit_encoding[i] : :(n)
        fun = sn == :(n) && multi_encodings ? :(alphabet{n}) : :(alphabet)
        argument = multi_encodings ? :(::Type{$type_name{$sn}}) : just_type
        push!(code.args, esc(:(function $fun($argument)
                               $(charset[i])
                               end)))
    end
    return code
end


@alphabet(DNAAlphabet,
          [2, 4],
          DNANucleotide,
          [DNA_A:DNA_T, alphabet(DNANucleotide)])


@alphabet(RNAAlphabet,
          [2, 4],
          RNANucleotide,
          [RNA_A:RNA_U, alphabet(RNANucleotide)])


@alphabet(AminoAcidAlphabet,
          8,
          AminoAcid,
          alphabet(AminoAcid))


# Encoders & Decoders
# -------------------

"""
Encode biological characters to binary representation.
"""
function encode end

"""
Decode biological characters from binary representation.
"""
function decode end

for (A, N) in ((DNAAlphabet, DNANucleotide),
               (RNAAlphabet, RNANucleotide)), n in (2, 4)
    @eval begin
        encode(::Type{$A{$n}}, x::$N) = reinterpret(UInt8, x)
        decode(::Type{$A{$n}}, x::UInt8) = reinterpret($N, x)
    end
end

encode(::Type{AminoAcidAlphabet}, x::AminoAcid) = reinterpret(UInt8, x)
decode(::Type{AminoAcidAlphabet}, x::UInt8) = reinterpret(AminoAcid, x)
