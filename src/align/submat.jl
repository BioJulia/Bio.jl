# Substitution Matrices
# ---------------------

"""
Supertype of substitution matrix.

The required method:

* `Base.getindex(submat, x, y) = <substitution score/cost from x to y>`
"""
abstract AbstractSubstitutionMatrix{T<:Real}


"""
Substitution matrix.
"""
type SubstitutionMatrix{T} <: AbstractSubstitutionMatrix{T}
    data::Matrix{T}
end

Base.convert(::Type{Matrix}, submat::SubstitutionMatrix) = submat.data
Base.convert{T}(::Type{Matrix{T}}, submat::SubstitutionMatrix{T}) = submat.data

@inline function Base.getindex(submat::SubstitutionMatrix, x, y)
    return submat.data[convert(UInt8, x)+1,convert(UInt8, y)+1]
end
function Base.setindex!(submat::SubstitutionMatrix, v, x, y)
    submat.data[convert(UInt8, x)+1,convert(UInt8, y)+1] = v
    return submat
end
function Base.fill!(submat::SubstitutionMatrix, v)
    return fill!(submat.data, v)
end


"""
Dichotomous substitution matrix.
"""
immutable DichotomousSubstitutionMatrix{T} <: AbstractSubstitutionMatrix{T}
    matching_score::T
    mismatching_score::T
end

@inline function Base.getindex(submat::DichotomousSubstitutionMatrix, x, y)
    return ifelse(
        x == y,
        submat.matching_score,
        submat.mismatching_score
    )
end

immutable UnitSubstitutionCost{T} <: AbstractSubstitutionMatrix{T} end
Base.getindex{T}(::UnitSubstitutionCost{T}, x, y) = ifelse(x == y, T(0), T(1))


# Utils for loading substitution matrices from files

function load_submat(name)
    submatfile = Pkg.dir("Bio", "src", "align", "data", "submat", name)
    return parse_ncbi_submat(submatfile)
end

function parse_ncbi_submat(filepath)
    ncols = 24
    submat = SubstitutionMatrix(Array{Int}(ncols, ncols))
    fill!(submat, 0)
    open(filepath) do io
        # column => one-letter amino acid
        header = Char[]
        # column width (inferred from the header line)
        w = 0
        for line in eachline(io)
            if startswith(line, "#")
                continue
            elseif isempty(header)
                # header line
                w = div(length(line) - 2, ncols)
                for col in 1:ncols
                    char = line[w*col+1]
                    @assert isalpha(char) || char == '*'
                    push!(header, char)
                end
            else
                a = line[1]
                if a in ('B', 'Z', '*')
                    # these amino acid code are not supported by Bio.Seq
                    continue
                end
                for col in 1:ncols
                    b = header[col]
                    if b in ('B', 'Z', '*')
                        # these amino acid code are not supported by Bio.Seq
                        continue
                    end
                    score = parse(Int, line[w*col-1:w*col+1])
                    submat[convert(AminoAcid, a), convert(AminoAcid, b)] = score
                end
            end
        end
    end
    return submat
end

# predefined substitution matrices
const PAM30  = load_submat("PAM30")
const PAM70  = load_submat("PAM70")
const PAM250 = load_submat("PAM250")
const BLOSUM45 = load_submat("BLOSUM45")
const BLOSUM50 = load_submat("BLOSUM50")
const BLOSUM62 = load_submat("BLOSUM62")
const BLOSUM80 = load_submat("BLOSUM80")
const BLOSUM90 = load_submat("BLOSUM90")
