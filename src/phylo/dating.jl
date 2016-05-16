module Dating

using Distributions, Roots

export DatingEstimate, coaltime


immutable DatingEstimate
    five::Float64
    fifty::Float64
    ninetyfive::Float64
end


function binomzero(p0::Float64, N::Int, B::Int)
    f(p::Float64) = cdf(Binomial(N, p), B) - p0
    return fzero(f, 0, 1)
end


function coaltime(len::Int, nmut::Int, mu::Float64)
    five = ceil(binomzero(0.05, len, nmut) / (2 * mu))
    fifty = ceil(binomzero(0.5, len, nmut) / (2 * mu))
    ninetyfive = ceil(binomzero(0.95, len, nmut) / (2 * mu))
    return DatingEstimate(five, fifty, ninetyfive)
end


end
