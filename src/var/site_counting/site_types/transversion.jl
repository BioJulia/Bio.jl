# transversion.jl
# ===============
#
# Define transition sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Transversion` site describes a site where two aligned nucleotides are
definately mutated, and the type of mutation is a transversion mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transversion mutation: i.e. A<->C,
A<->T, G<->T, G<->C.
"""
immutable Transversion <: Mutation end

# Methods for the naive framework.
# --------------------------------

for Alph in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function ischange{A<:$Alph,B<:$Alph}(::Type{Transversion}, a::BioSequence{A}, b::BioSequence{B}, idx)
            ai = a[idx]
            bi = b[idx]
            return (ai != bi) & ((ispurine(ai) & ispyrimidine(bi)) | (ispyrimidine(ai) & ispurine(bi)))
        end
    end
end

# Methods for the bitparallel framework.
# --------------------------------------
