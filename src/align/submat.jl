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
    submat::Matrix{T}
end

@inline function Base.getindex(submat::SubstitutionMatrix, x, y)
    return submat.submat[convert(UInt8, x)+1,convert(UInt8, y)+1]
end
function Base.setindex!(submat::SubstitutionMatrix, v, x, y)
    submat.submat[convert(UInt8, x)+1,convert(UInt8, y)+1] = v
    return submat
end
function Base.fill!(submat::SubstitutionMatrix, v)
    return fill!(submat.submat, v)
end


"""
Dichotomous substitution matrix.
"""
immutable DichotomousSubstitutionMatrix{T} <: AbstractSubstitutionMatrix{T}
    matching_score::T
    mismatching_score::T
end

function Base.getindex(submat::DichotomousSubstitutionMatrix, x, y)
    return ifelse(
        x == y,
        submat.matching_score,
        submat.mismatching_score
    )
end

function load_submat(name)
    submatfile = Pkg.dir("Bio", "src", "align", "data", "submat", name)
    return parse_ncbi_submat(submatfile)
end

function parse_ncbi_submat(filepath)
    n = 24
    submat = SubstitutionMatrix(Array{Int}(n, n))
    fill!(submat, 0)
    open(filepath) do io
        # column => one-letter amino acid
        header = Char[]
        for line in eachline(io)
            if startswith(line, "#")
                continue
            elseif isempty(header)
                for col in 1:n
                    char = line[3col+1]
                    @assert char != ' '
                    push!(header, char)
                end
            else
                a = line[1]
                if a in ('B', 'Z', '*')
                    # these amino acid code are not supported by Bio.Seq
                    continue
                end
                for col in 1:n
                    b = header[col]
                    if b in ('B', 'Z', '*')
                        # these amino acid code are not supported by Bio.Seq
                        continue
                    end
                    score = parse(Int, line[3col-1:3col+1])
                    submat[convert(AminoAcid, a), convert(AminoAcid, b)] = score
                end
            end
        end
    end
    return submat
end

const BLOSUM62 = load_submat("BLOSUM62")
