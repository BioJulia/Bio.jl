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
       middle, lower

"""
An abstract dating estimate type.

Types inheriting from DatingEstimate should have the following methods defined:

* lower: Returns the lower bound of the dating estimate.
* upper: Returns the upper bound of the dating estimate.

They may define their own show methods.
"""
abstract DatingEstimate
abstract DatingMethod

Base.print(io::IO, de::DatingEstimate) = println(io, "$(lower(de)) ... $(upper(de))")
Base.show(io::IO, de::DatingEstimate) = println(io, "Coalesence time estimate between two sequences:\nt lies between $(lower(de)) and $(upper(de))")
Base.in(val::Float64, de::DatingEstimate) = val <= upper(de) && val >= lower(de)

"""
A very simple expected divergence time estimate.
Assumes a strict molecult clock and that the divergence time is equal to

\$ t = d / (2\\mu) \$

Where \$d\$ is the evolutionary distance computed for two aligned sequences,
and \$\\mu\$ is the substitution rate.
"""
immutable SimpleEstimate <: DatingMethod end
restype(::Type{SimpleEstimate}) = Float64

"""
    coaltime(N::Int, K::Int, µ::Float64, ::Type{SimpleEstimate})

Compute the coalescence time between two sequences by the `SimpleEstimate`
method.

`N` is the length of the two aligned sequences, `K` is the number of mutations
and `µ` is the assumed mutation rate.

# Examples
```julia
coaltime(50, 17, 10e-9, SimpleEstimate)
```
"""
@inline function coaltime(N::Int, K::Int, µ::Float64, ::Type{SimpleEstimate})
    @assert N >= K >= 0 error("Condition: N >= K >= 0, not met.")
    return (K / N) / (2 * µ)
end

"""
    coaltime(N::Int, d::Float64, µ::Float64, ::Type{SimpleEstimate})

Compute the coalescence time between two sequences by the `SimpleEstimate`
method.

`N` is the length of the two aligned sequences, `d` is the evolutionary distance
and `µ` is the assumed mutation rate.

# Examples
```julia
coaltime(50, 0.12, 10e-9, SimpleEstimate)
```
"""
@inline function coaltime(N::Int, d::Float64, µ::Float64, ::Type{SimpleEstimate})
    @assert N >= d >= 0 error("Genetic distance `d` must be a value between 1 and 0.")
    return d / (2 * µ)
end

"""
`SpeedDate` is the name given to a method of estimating a divergence time between
two DNA sequence regions that was first implemented in the R package
HybridCheck in order to date regions of introgression in large sequence
contigs.

The coalescence time is estimated using the number of mutations that have
occurred between two aligned sequences. The calculation uses a strict molecular
clock which assumes a constant substitution rate, both through time and across
taxa.
Modelling the mutation accumulation process as a Bernoulli trial, the
probability of observing  \$k\$  or fewer mutations between two sequences of
length \$n\$ can be given as:

\$ Pr(X \\le k) = \\sum_{i=0}^{\\lfloor k \\rfloor} \\binom{n}{i} p^i (1 - p)^{n-i} \$

Where \$p\$ is the probability of observing a single mutation between the two
aligned sequences.
The value of \$p\$ depends on two key factors: the substitution rate and the
coalescence time.
If you assume a molecular clock, whereby two DNA sequences are both accumulating
mutations at a rate \$\\mu\$ for \$t\$ generations, then you may define
\$p = 2\\mu t\$.

Using these assumptions, the SpeedDate method finds the root of the following
formula for \$Pr(X \\le k) = 0.05\$, \$0.5\$, and \$0.95\$, and then divides the
three answers by twice the assumed substitution rate.

\$ f(n, k, 2\\mu t, Pr(X \\le k) = \\left( \\sum_{i=0}^{\\lfloor k \\rfloor} \\binom{n}{i} {2\\mu t}^i (1 - 2\\mu t)^{n-i}   \\right) - Pr(X \\le k) \$

This results in an upper, middle, and lower estimate of the coalescence time
\$t\$ of the two sequences (expressed as the number of generations).
"""
immutable SpeedDating <: DatingMethod end
restype(::Type{SpeedDating}) = SDResult

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
middle(x::Float64) = x
middle(x::SDResult) = x.middle
upper(c::Float64) = x
upper(x::SDResult) = x.upper

function Base.show(io::IO, de::SDResult)
    println(io, "Coalescence time estimate:\n5%: $(de.lower), 95%: $(de.upper)")
end

@inline function binomzero(p0::Float64, N::Int, B::Int)
    f(p::Float64) = cdf(Binomial(N, p), B) - p0
    return fzero(f, 0, 1)
end

"""
    coaltime(N::Int, K::Int, µ::Float64, ::Type{SpeedDating})

Compute the coalescence time between two sequences by modelling the process of
mutation accumulation between two sequences as a bernoulli process.

This method was first described in the paper:
Ward, B. J., & van Oosterhout, C. (2016). Hybridcheck: Software for the rapid detection, visualization and dating of recombinant regions in genome sequence data. Molecular Ecology Resources, 16(2), 534–539.

`N` is the length of the two aligned sequences, `K` is the number of mutations
and `µ` is the assumed mutation rate.

# Examples
```julia
coaltime(50, 17, 10e-9, SpeedDating)
```
"""
@inline function coaltime(N::Int, K::Int, µ::Float64, ::Type{SpeedDating})
    @assert N >= K >= 0 error("Condition: N >= K >= 0, not met.")
    div = 2 * µ
    ninetyfive = ceil(binomzero(0.05, N, K) / div)
    fifty = ceil(binomzero(0.5, N, K) / div)
    five = ceil(binomzero(0.95, N, K) / div)
    return SDResult(five, fifty, ninetyfive)
end

"""
    coaltime(N::Int, p::Float64, µ::Float64, ::Type{SpeedDating})

Compute the coalescence time between two sequences by modelling the process of
mutation accumulation between two sequences as a bernoulli process.

This method was first described in the paper:
Ward, B. J., & van Oosterhout, C. (2016). Hybridcheck: Software for the rapid
detection, visualization and dating of recombinant regions in genome sequence
data. Molecular Ecology Resources, 16(2), 534–539.

`N` is the length of the two aligned sequences, `d` is the evolutionary distance
and `µ` is the assumed mutation rate.

# Examples
```julia
coaltime(50, 17, 10e-9, SpeedDating)
```
"""
@inline function coaltime(N::Int, d::Float64, µ::Float64, ::Type{SpeedDating})
    @assert 1 >= d >= 0 error("Genetic distance `d` must be a value between 1 and 0.")
    K = convert(Int, ceil(d * N))
    return coaltime(N, K, µ, SpeedDating)
end

"""
    coaltime{M<:DatingMethod,R<:Real}(N::Int, arr::AbstractArray{R}, µ::Float64, ::Type{M})

Compute the coalescence time between two sequences by modelling the process of
mutation accumulation between two sequences as a bernoulli process.

This method was first described in the paper:
Ward, B. J., & van Oosterhout, C. (2016). Hybridcheck: Software for the rapid detection, visualization and dating of recombinant regions in genome sequence data. Molecular Ecology Resources, 16(2), 534–539.

In this specific method, `N` is the length of the aligned sequences,
`arr` is an array of either mutation counts (Integers) or genetic distances
(Floats), `µ` is the assumed mutation rate.

# Examples
```julia
coaltime(50, [17, 20, 10, 7], 10e-9, SpeedDating)
coaltime(50, [0.01, 0.09, 0.12, 0.20], 10e-9, SpeedDating)
```
"""
@inline function coaltime{M<:DatingMethod,R<:Real}(N::Int,
                                           arr::AbstractArray{R},
                                           µ::Float64,
                                           ::Type{M})
    Ts = similar(arr, restype(M))
    @inbounds for i in eachindex(arr)
        Ts[i] = coaltime(N, arr[i], µ, M)
    end
    return Ts
end

end
