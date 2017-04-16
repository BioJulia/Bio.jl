# site_types.jl
# =============
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
An abstract type which you can inherit from to build on the site counting
"""
abstract Site

@inline function count_algorithm(s::Site, a::BioSequence, b::BioSequence)
    return NAIVE
end

# For a given site type you must define the following functions for the
# generated function and method dispatch to work on:

# Methods required for both the naive and bitparallel algorithms.
# ---------------------------------------------------------------

# What will the output type of count(YOUR_SITE_TYPE, seqa, seqb), be?
# Note that the type specified here, should have a Base.zero or start_counter
# method defined for it.
counter_type{S<:Site,A<:Alphabet}(::Type{S}, ::Type{A}) = Int
counter_type{S<:Site,A<:Alphabet,B<:Alphabet}(::Type{S}, ::Type{A}, ::Type{B}) = Int

# How to start the count accumulator.
# The default method is to use Base.zero on the `counter_type`, but you can overload
# this.
start_counter{S<:Site,A<:Alphabet}(::Type{S}, ::Type{A}) = zero(counter_type(S, A))
start_counter{S<:Site,A<:Alphabet,B<:Alphabet}(::Type{S}, ::Type{A}, ::Type{B}) = zero(counter_type(S, A, B))

# Methods required for just the naive algorithm.
# ----------------------------------------------

# Methods to update the count when given the output of the `issite` method you
# define for your site type.
# The naieve algorithm loops over every site and applies the `issite` method,
# typically yielding either `true` or `false`. But by defining your own
# `start_count`, `update_count`, and `issite` methods this can be different.
update_counter(acc::Int, up::Bool) = acc + up

# Methods required for the bitparallel algorithm.
# -----------------------------------------------

# Methods to update the count when given the output of the `count_bitpairs`
# and `count_nibbles` methods you define for your site type.
# The bitparallel algorithm loops over whole words of binary data,
# and applies the `count_bitpairs` method, typically yielding and integer.
# But again, by defining your own `start_counter`, `update_counter` and
# `count_bitpairs`
update_counter(acc::Int, up::Int) = acc + up

include("certain.jl")
include("ambiguous.jl")
include("gap.jl")
include("match.jl")
include("mismatch.jl")
