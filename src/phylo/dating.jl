module Dating

using Distributions, Roots

export DatingEstimate, coaltime

"""
DatingEstimate is a simple datatype to store the results of the coaltime
function.

coaltime estimates date range that may be considered a 95% confidence range in
which the true coalsescence time between two sequences lies.

The type store the 5%, 50%, and 95% values for the range.
"""
immutable DatingEstimate
    five::Float64
    fifty::Float64
    ninetyfive::Float64
end


function binomzero(p0::Float64, N::Int, B::Int)
    f(p::Float64) = cdf(Binomial(N, p), B) - p0
    return fzero(f, 0, 1)
end

"""
    coaltime(len::Int, nmut::Int, mu::Float64)

Compute the coalescence time between two sequences by modelling the process of
mutation accumulation between two sequences as a bernoulli process.

This method was first described in the paper:
Ward, B. J., & van Oosterhout, C. (2016). Hybridcheck: Software for the rapid detection, visualization and dating of recombinant regions in genome sequence data. Molecular Ecology Resources, 16(2), 534–539.

len is the length of the two aligned sequences, nmut is the number of mutations
and mu is the assumed mutation rate.

# Examples
```julia
coaltime(50, 17, 10e-9)
```
"""
function coaltime(len::Int, nmut::Int, mu::Float64)
    five = ceil(binomzero(0.05, len, nmut) / (2 * mu))
    fifty = ceil(binomzero(0.5, len, nmut) / (2 * mu))
    ninetyfive = ceil(binomzero(0.95, len, nmut) / (2 * mu))
    return DatingEstimate(five, fifty, ninetyfive)
end


end
