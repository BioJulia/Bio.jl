# phylo/dating.jl
# ==================
#
# Types and methods for computing coalescence times between two sequences.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md



module Dating

using Distributions, Roots

export DatingEstimate, coaltime, SimpleEstimate, SpeedDating, SDResult, upper,
       lower

"""
An abstract dating estimate type.

Types inheriting from DatingEstimate should have the following methods defined:

* lower: Returns the lower bound of the dating estimate.
* upper: Returns the upper bound of the dating estimate.

They may define their own show methods.
"""
abstract DatingEstimate

Base.print(io::IO, de::DatingEstimate) = println(io, "$(lower(de)) ... $(upper(de))")
Base.show(io::IO, de::DatingEstimate) = println(io, "Coalesence time estimate between two sequences:\nt lies between $(lower(de)) and $(upper(de))")
Base.in(val::Float64, de::DatingEstimate) = val <= upper(de) && val >= lower(de)

# Different coalescence time estimation algorithms.
# method coaltime is dispatched according to arguments or type of
immutable SimpleEstimate end

"""
    coaltime(len::Int, nmut::Int, mu::Float64, ::Type{SimpleEstimate})

Compute the coalescence time between two sequences by assuming a mutation rate,
and then assuming the divergence between two sequences is 2 * mu * t

# Examples
```julia
coaltime(50, 17, 10e-9, SimpleEstimate)
```
"""
function coaltime(len::Int, nmut::Int, mu::Float64, ::Type{SimpleEstimate})
    return (nmut / len) / (2 * mu)
end


immutable SpeedDating end

"""
SDResult is a simple datatype to store the results of the coaltime
function.

coaltime estimates date range that may be considered a 95% confidence range in
which the true coalsescence time between two sequences lies.

The type store the 5%, 50%, and 95% values for the range.
"""
immutable SDResult <: DatingEstimate
    lower::Float64
    middle::Float64
    upper::Float64
end

lower(x::Float64) = x
lower(x::SDResult) = x.lower
upper(c::Float64) = x
upper(x::SDResult) = x.upper

function Base.show(io::IO, de::SDResult)
    println(io, "Coalescence time estimate:\n5%: $(de.lower), 95%: $(de.upper)")
end

function binomzero(p0::Float64, N::Int, B::Int)
    f(p::Float64) = cdf(Binomial(N, p), B) - p0
    return fzero(f, 0, 1)
end

"""
    coaltime(len::Int, nmut::Int, mu::Float64, ::Type{SpeedDating})

Compute the coalescence time between two sequences by modelling the process of
mutation accumulation between two sequences as a bernoulli process.

This method was first described in the paper:
Ward, B. J., & van Oosterhout, C. (2016). Hybridcheck: Software for the rapid detection, visualization and dating of recombinant regions in genome sequence data. Molecular Ecology Resources, 16(2), 534–539.

len is the length of the two aligned sequences, nmut is the number of mutations
and mu is the assumed mutation rate.

# Examples
```julia
coaltime(50, 17, 10e-9, SpeedDating)
```
"""
function coaltime(len::Int, nmut::Int, mu::Float64, ::Type{SpeedDating})
    ninetyfive = ceil(binomzero(0.05, len, nmut) / (2 * mu))
    fifty = ceil(binomzero(0.5, len, nmut) / (2 * mu))
    five = ceil(binomzero(0.95, len, nmut) / (2 * mu))
    return SDResult(five, fifty, ninetyfive)
end


end
