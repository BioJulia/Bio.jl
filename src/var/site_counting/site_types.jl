# Different types of site that can be counted
# ===========================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract Site
abstract Mutation <: Site

"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol or a gap.
"""
immutable Certain <: Site end

"""
An `Ambiguous` site describes a site where either of two aligned sites are an
ambiguity symbol.
"""
immutable Ambiguous <: Site end

"""
An `Gap` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
immutable Gap <: Site end

"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
immutable Match <: Site end

"""
A `Mismatch` site describes a site where two aligned nucleotides are not the
same biological symbol.
"""
immutable Mismatch <: Site end

"""
A `Mismatch` site describes a site where two aligned nucleotides are definately
conserved. By definately conserved this means that the symbols of the site are
non-ambiguity symbols, and they are the same symbol.
"""
immutable Conserved <: Mutation end

"""
A `Mutated` site describes a site where two aligned nucleotides are definately
mutated. By definately mutated this means that the symbols of the site are
non-ambiguity symbols, and they are not the same symbol.
"""
immutable Mutated <: Mutation end

"""
A `Transition` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transition mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transition mutation: i.e. A<->G, or C<->T.
"""
immutable Transition <: Mutation end

"""
A `Transversion` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transversion mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transversion mutation: i.e. A<->C,
A<->T, G<->T, G<->C.
"""
immutable Transversion <: Mutation end
