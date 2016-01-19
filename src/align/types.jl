# Alignment Types
# ---------------

abstract AbstractAlignment

"""Global-global alignment with end gap penalties."""
immutable GlobalAlignment <: AbstractAlignment end

"""Global-local alignment."""
immutable SemiGlobalAlignment <: AbstractAlignment end

"""Global-global alignment without end gap penalties."""
immutable OverlapAlignment <: AbstractAlignment end

"""Local-local alignment."""
immutable LocalAlignment <: AbstractAlignment end


# distances

"""Edit distance."""
immutable EditDistance <: AbstractAlignment end

"""
Levenshtein distance.

A special case of `EditDistance` with the costs of mismatch, insertion, and
deletion are 1.
"""
immutable LevenshteinDistance <: AbstractAlignment end

"""
Hamming distance.

A special case of `EditDistance` with the costs of insertion and deletion are
infinitely large.
"""
immutable HammingDistance <: AbstractAlignment end
