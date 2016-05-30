
#=
Objects of type BinarySNPs store the SNPs of an individual
The name BinarySNPs is double edged: it refers to both the fact that
SNPs are stored in a binary format, and also that the SNPs are binary,
in that they have two states, a reference state and an alternate state.
=#

type BinarySNPs
    snps::Vector{Vector{UInt8}}
    numberOfLoci::Int
    naLoci::Vector{Int}
    label::ASCIIString
    ploidy::Int
end

function BinarySNPs(snps::DataArray{Int}, label::ASCIIString = "")
    nas, snps = resolvearray(snps)
    return BinarySNPs(snps, nas, label)
end

function BinarySNPs(snps::DataArray{Int}, ploidy::Int, label::ASCIIString = "")
    nas, snps = resolvearray(snps)
    return BinarySNPs(snps, nas, ploidy, label)
end

# Automatically determine the ploidy from the input vector.
function BinarySNPs(snps::Vector{Int}, na::Vector{Int}, label::ASCIIString = "")
    ploidy = maximum(snps)
    if ploidy <= 0
        ploidy = 1
    end
    return BinarySNPs(snps, na, ploidy, label)
end

function BinarySNPs(snps::Vector{Int}, na::Vector{Int}, ploidy::Int, label::ASCIIString = "")
    snps[na] = 0
    @assert length(snps) > 0
    @assert ploidy > 0
    for snp in snps
        if snp < 0 || snp > ploidy
            error("SNP values cannot be less than 0 i.e. cannot be negative.")
        end
    end
    maxallele = maximum(snps)
    byteSNPs = Vector{Vector{UInt8}}(maxallele)
    if maxallele > 1
        snpListInd = 0
        while maxallele > 0
            snpListInd += 1
            temporary = Vector{Int}(snps .== maxallele)
            byteSNPs[snpListInd] = binary_integers_to_raw(temporary)
            snps = snps - temporary
            maxallele = maximum(snps)
        end
    else
        byteSNPs[1] = binary_integers_to_raw(snps)
    end
    return BinarySNPs(byteSNPs, length(snps), na, label, ploidy)
end

"""
Internal function which converts vectors of integers storing 1's and 0's, into
a vector of bytes which stores each number as a single bit.
"""
function binary_integers_to_raw(snps::Vector{Int})
    for snp in snps
        if snp > 1 || snp < 0
            error("Input vector contained numbers other than 1 or 0.")
        end
    end
    nbytes = div(length(snps), 8)
    nbytes += mod(length(snps), 8) > 0 ? 1 : 0
    newlength = 8 * nbytes
    snps = vcat(snps, zeros(Int, newlength - length(snps)))
    bytevector = Vector{UInt8}(nbytes)
    fill!(bytevector, 0)
    byteidx = 1
    bitidx = 1
    for snpidx in 1:length(snps)
        bytevector[byteidx] += UInt8((2 ^ (bitidx - 1)) * snps[snpidx])
        if bitidx == 8
            byteidx += 1
            bitidx = 1
        else
            bitidx += 1
        end
    end
    return bytevector
end

"""
Internal functions for writing the bits of a byte, to a vector of integers.
Will not modify snps, but will modify out.
"""
function raw_to_binary_integers!(byte::UInt8, vec::Vector{Int})
    rest = Int(byte)
    fill!(vec, 0)
    for idx in 8:-1:1
        ref = 2 ^ (idx - 1)
        if rest >= ref
            vec[idx] = 1
            rest -= ref
            if rest == 0
                break
            end
        end
    end
end

function raw_to_binary_integers(bytes::Vector{UInt8})
    temporary = Vector{Int}(8)
    result = Vector{Int}(length(bytes) * 8)
    residx = 0
    for byteidx in 1:length(bytes)
        raw_to_binary_integers!(bytes[byteidx], temporary)
        for j in 1:8
            result[j + residx] = temporary[j]
        end
        residx += 8
    end
    return result
end

"Internal function used to split DataArrays{Int} into two vectors."
function resolvearray(snps::DataArray{Int})
    return find(isna(snps)), snps.data
end

nacount(snps::BinarySNPs) = length(snps.naLoci)

function Base.show(io::IO, snps::BinarySNPs)
    if snps.label != ""
        write(io, "$(snps.label)\n")
    end
    write(io, "$(snps.numberOfLoci) Binary SNPs stored in compact bit format.\n")
    write(io, "Ploidy: $(snps.ploidy)\n")
    numna = nacount(snps)
    write(io, "$(numna) loci are NA. ($(round((numna / snps.numberOfLoci) * 100))%)\n")
    write(io, "Size: $((length(snps.snps[1]) * length(snps.snps)) / 1000) Kb\n")
end

function Base.convert(::Type{Vector{Int}}, binsnps::BinarySNPs)
    snps = binsnps.snps
    nbytes = length(snps[1])
    combnvec = vcat(snps...)
    result = Vector{Int}(nbytes * 8)
    fill!(result, 0)
    temporary = Vector{Int}(8)
    for vecidx in 0:length(snps) - 1
        idres = 0
        for byteidx in 0:nbytes - 1
            raw_to_binary_integers!(combnvec[byteidx + (vecidx * nbytes) + 1], temporary)
            for tempval in 1:8
                result[tempval + idres] += temporary[tempval]
            end
            idres += 8
        end
    end
    return result[1:binsnps.numberOfLoci]
end

function Base.convert(::Type{DataArray{Int}}, binsnps::BinarySNPs)
    arr = DataArray{Int}(Vector{Int}(binsnps))
    arr[binsnps.naLoci] = NA
    return arr
end
