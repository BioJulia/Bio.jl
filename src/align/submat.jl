# Substitution Matrices
# =====================
#
# Substitution matrix types and predefined matrices.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
Supertype of substitution matrix.

The required method:

* `Base.getindex(submat, x, y)`: substitution score/cost from `x` to `y`
"""
abstract AbstractSubstitutionMatrix{S<:Real}

"""
Substitution matrix.
"""
immutable SubstitutionMatrix{T,S} <: AbstractSubstitutionMatrix{S}
    # square substitution matrix
    data::Matrix{S}
    # score is defined or not
    defined::BitMatrix

    function SubstitutionMatrix(data::Matrix{S}, defined::BitMatrix)
        @assert size(data) == size(defined)
        @assert size(data, 1) == size(data, 2) == length(Bio.Seq.alphabet(T)) - 1
        return new(data, defined)
    end

end

function SubstitutionMatrix{T,S}(::Type{T},
                                 submat::AbstractMatrix{S},
                                 defined::AbstractMatrix{Bool},
                                 default_match::S,
                                 default_mismatch::S)
    alpha = alphabet_without_gap(T)
    n = length(alpha)
    data = Matrix{S}(n, n)
    for x in alpha, y in alpha
        i = index(x)
        j = index(y)
        if defined[i,j]
            data[i,j] = submat[i,j]
        else
            data[i,j] = x == y ? default_match : default_mismatch
        end
    end
    return SubstitutionMatrix{T,S}(data, copy(defined))
end

function SubstitutionMatrix{T,S}(::Type{T},
                                 submat::AbstractMatrix{S},
                                 # assume scores of all substitutions are defined
                                 defined=trues(size(submat));
                                 default_match=S(0), default_mismatch=S(0))
    return SubstitutionMatrix(T, submat, defined, default_match, default_mismatch)
end

function SubstitutionMatrix{T,S}(scores::Associative{Tuple{T,T},S};
                                 default_match=S(0), default_mismatch=S(0))
    n = length(Bio.Seq.alphabet(T)) - 1
    submat = Matrix{S}(n, n)
    defined = falses(n, n)
    for ((x, y), score) in scores
        i = index(x)
        j = index(y)
        submat[i,j] = score
        defined[i,j] = true
    end
    return SubstitutionMatrix(T, submat, defined, default_match, default_mismatch)
end

Base.convert(::Type{Matrix}, submat::SubstitutionMatrix) = copy(submat.data)

@inline function Base.getindex{T}(submat::SubstitutionMatrix{T}, x, y)
    i = index(convert(T, x))
    j = index(convert(T, y))
    return submat.data[i,j]
end

function Base.setindex!{T}(submat::SubstitutionMatrix{T}, val, x, y)
    i = index(convert(T, x))
    j = index(convert(T, y))
    submat.data[i,j] = val
    submat.defined[i,j] = true
    return submat
end

function Base.copy{T,S}(submat::SubstitutionMatrix{T,S})
    return SubstitutionMatrix{T,S}(copy(submat.data), copy(submat.defined))
end

Base.minimum(submat::SubstitutionMatrix) = minimum(submat.data)
Base.maximum(submat::SubstitutionMatrix) = maximum(submat.data)

function Base.show{T,S}(io::IO, submat::SubstitutionMatrix{T,S})
    alpha = alphabet_without_gap(T)
    n = length(alpha)
    mat = Matrix{String}(n, n)
    for (i, x) in enumerate(alpha), (j, y) in enumerate(alpha)
        i′ = index(x)
        j′ = index(y)
        mat[i,j] = string(submat.data[i′,j′])
        if !submat.defined[i′,j′]
            mat[i,j] = underline(mat[i,j])
        end
    end

    # add rows and columns
    rows = map(string, alpha)
    cols = vcat("", rows)
    mat = hcat(rows, mat)
    mat = vcat(reshape(cols, 1, length(cols)), mat)

    println(io, summary(submat), ':')
    width = maximum(map(x -> length(x), mat))
    for i in 1:n+1
        for j in 1:n+1
            print(io, lpad(mat[i,j], width + 1))
        end
        println(io)
    end
    print(io, "(underlined values are default ones)")
end

underline(s) = join([string(c, '\U0332') for c in s])

# Return a vector of all symbols of `T` except the gap symbol.
function alphabet_without_gap{T}(::Type{T})
    return filter!(x -> x != Bio.Seq.gap(T), collect(Bio.Seq.alphabet(T)))
end

# Return the row/column index of `nt`.
function index(nt::Bio.Seq.NucleicAcid)
    return convert(Int, nt)
end

# Return the row/column index of `aa`.
function index(aa::Bio.Seq.AminoAcid)
    return convert(Int, aa) + 1
end


"""
Dichotomous substitution matrix.
"""
immutable DichotomousSubstitutionMatrix{S} <: AbstractSubstitutionMatrix{S}
    match::S
    mismatch::S
end

function DichotomousSubstitutionMatrix(match::Real, mismatch::Real)
    typ = promote_type(typeof(match), typeof(mismatch))
    return DichotomousSubstitutionMatrix{typ}(match, mismatch)
end

function Base.convert{T,S}(::Type{SubstitutionMatrix{T,S}},
                           submat::DichotomousSubstitutionMatrix)
    n = length(Bio.Seq.alphabet(T)) - 1
    data = Matrix{S}(n, n)
    fill!(data, submat.mismatch)
    data[diagind(data)] = submat.match
    return SubstitutionMatrix{T,S}(data, trues(n, n))
end

@inline function Base.getindex(submat::DichotomousSubstitutionMatrix, x, y)
    return ifelse(x == y, submat.match, submat.mismatch)
end

function Base.show(io::IO, submat::DichotomousSubstitutionMatrix)
    match = string(submat.match)
    mismatch = string(submat.mismatch)
    w = max(length(match), length(mismatch))
    println(io, typeof(submat), ':')
    # right aligned
    print(io, "     match = ", lpad(match, w), '\n')
    print(io, "  mismatch = ", lpad(mismatch, w))
end


# Predefined Substitution Matrices
# --------------------------------

function load_submat{T}(::Type{T}, name)
    return parse_ncbi_submat(T, joinpath(dirname(@__FILE__), "data", "submat", name))
end

function parse_ncbi_submat{T}(::Type{T}, filepath)
    scores = Dict{Tuple{T,T},Int}()
    cols = T[]
    open(filepath) do io
        for line in eachline(io)
            line = chomp(line)
            if startswith(line, "#") || isempty(line)
                continue  # skip comments and empty lines
            end
            if isempty(cols)
                for y in matchall(r"[A-Z*]", line)
                    @assert length(y) == 1
                    push!(cols, convert(T, first(y)))
                end
            else
                x = convert(T, first(line))
                ss = matchall(r"-?\d+", line)
                @assert length(ss) == length(cols)
                for (y, s) in zip(cols, ss)
                    scores[(x, y)] = parse(Int, s)
                end
            end
        end
    end
    return SubstitutionMatrix(scores, default_match=0, default_mismatch=0)
end

const EDNAFULL = load_submat(Bio.Seq.DNA, "NUC.4.4")
const PAM30    = load_submat(Bio.Seq.AminoAcid, "PAM30")
const PAM70    = load_submat(Bio.Seq.AminoAcid, "PAM70")
const PAM250   = load_submat(Bio.Seq.AminoAcid, "PAM250")
const BLOSUM45 = load_submat(Bio.Seq.AminoAcid, "BLOSUM45")
const BLOSUM50 = load_submat(Bio.Seq.AminoAcid, "BLOSUM50")
const BLOSUM62 = load_submat(Bio.Seq.AminoAcid, "BLOSUM62")
const BLOSUM80 = load_submat(Bio.Seq.AminoAcid, "BLOSUM80")
const BLOSUM90 = load_submat(Bio.Seq.AminoAcid, "BLOSUM90")
