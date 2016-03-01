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
macro alphabet(args...)
    code = quote end

    # Process options for desired alphabet
    # -------------------------------------

    options = Dict{Symbol, Any}()
    for arg in args
        options[arg.args[1]] = arg.args[2]
    end

    # Check for and process base name for alphabet type.
    if !haskey(options, :name)
        error("No name provided for alphabet.")
    end
    typename = options[:name]

    # Check for and process whether the alphabet will have
    # different bit encodings.
    bits = get(options, :bits, 8)
    bits = typeof(bits) == Expr ? bits.args : bits

    # Check for and process different element types.
    if !haskey(options, :element_type)
        error("You need to provide an element type.")
    end
    element_type = options[:element_type]
    element_type = typeof(element_type) == Expr ? element_type.args : element_type


    # Check for and make sure that the alphabet argument has
    # been provided.
    if !haskey(options, :alphabet)
        error("You must provide one or more alphabets.")
    end
    alphs = options[:alphabet]
    alphs = typeof(alphs) == Expr ? alphs.args : alphs


    # Make some descision on how functions and code need to be defined
    # -----------------------------------------------------------------

    # Descisions that depend on if alphabet has different bit encodings:
    typename_n = esc(:($typename))
    put_n = false
    if typeof(bits) <: Array
        if length(bits) > 1
            typename_n = esc(:($typename{n}))
            put_n = true
        else
            bits = bits[1]
        end
    end

    # Descisions that depend on arguments passed as alphabet:
    multi_eltypes = false
    if typeof(element_type) <: Array
        multi_eltypes = length(element_type) > 1
        if multi_eltypes && length(element_type) != length(bits)
            error("Number of element options does not match the number of bit options.")
        end
    else
        element_type = Symbol[element_type]
    end


    # Descisions based on arguments passed as alphabet:
    multi_alphs = false
    if typeof(alphs) <: Array
        # Determine if it's one array of a single type, an array of character sets.
        types = Vector{DataType}(length(alphs))
        for i in 1:length(alphs)
            if typeof(alphs[i]) == Expr
                if alphs[i].head == :(:)
                    types[i] = Range
                elseif alphs[i].head == :(vect)
                    types[i] = Array
                elseif alphs[i].head == :(call)
                    types[i] = Function
                else
                    error("Unrecognised head when determining types of alphabet arguments.")
                end
            else
                types[i] = typeof(alphs[i])
            end
        end
        multi_alphs = all(Bool[i == Range || i == Array || Function for i in types]) && length(types) > 1
        if multi_alphs
            if length(alphs) != length(bits)
                error("Number of alphabet options does not match number of bit options.")
            end
        else
            alphs = Array{Int, 1}[[alphs]]
        end
    end


    # Build up code expressions now
    # ------------------------------

    # Type definition:
    push!(code.args, :(immutable $typename_n <: Alphabet end))

    # bitsof functions:
    for i in bits
        t = put_n ? :(::Type{$typename{$i}}) : :(::Type{$typename})
        push!(code.args, esc(:(function bitsof($t)
                               $i
                               end)))
    end

    # Base.eltype functions:
    for i in 1:length(element_type)
        n = multi_eltypes ? bits[i] : :(n)
        t = put_n ? :(::Type{$typename{$n}}) : :(::Type{$typename})
        push!(code.args, esc(:(function Base.eltype($t)
                               $(element_type[i])
                               end)))
    end

    # alphabet functions:
    for i in 1:length(alphs)
        n = multi_alphs ? bits[i] : :(n)
        t = put_n ? :(::Type{$typename{$n}}) : :(::Type{$typename})
        push!(code.args, esc(:(function alphabet($t)
                               $(alphs[i])
                               end)))
    end

    return code
end


"""
DNA nucleotide alphabet.
"""
#immutable DNAAlphabet{n} <: Alphabet end

@alphabet(name = DNAAlphabet,
          bits = [2, 4],
          element_type = DNANucleotide,
          alphabet = [DNA_A:DNA_T, alphabet(DNANucleotide)])

"""
RNA nucleotide alphabet.
"""
#immutable RNAAlphabet{n} <: Alphabet end

@alphabet(name = RNAAlphabet,
          bits = [2, 4],
          element_type = RNANucleotide,
          alphabet = [RNA_A:RNA_U, alphabet(RNANucleotide)])

"""
Amino acid alphabet.
"""
#immutable AminoAcidAlphabet <: Alphabet end
@alphabet(name = AminoAcidAlphabet,
          bits = 8,
          element_type = AminoAcid,
          alphabet = alphabet(AminoAcid))


#for n in (2, 4)
#    @eval begin
#        bitsof(::Type{DNAAlphabet{$n}}) = $n
#        bitsof(::Type{RNAAlphabet{$n}}) = $n
#    end
#end
#bitsof(::Type{AminoAcidAlphabet}) = 8

#Base.eltype(::Type{DNAAlphabet}) = DNANucleotide
#Base.eltype(::Type{RNAAlphabet}) = RNANucleotide
#Base.eltype{n}(::Type{DNAAlphabet{n}}) = DNANucleotide
#Base.eltype{n}(::Type{RNAAlphabet{n}}) = RNANucleotide
#Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid

#alphabet(::Type{DNAAlphabet{2}}) = DNA_A:DNA_T
#alphabet(::Type{RNAAlphabet{2}}) = RNA_A:RNA_U
#alphabet(::Type{DNAAlphabet{4}}) = alphabet(DNANucleotide)
#alphabet(::Type{RNAAlphabet{4}}) = alphabet(RNANucleotide)
#alphabet(::Type{AminoAcidAlphabet}) = alphabet(AminoAcid)


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
