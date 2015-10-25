# Substitution Matrices
# ---------------------

"""
Supertype of substitution matrix.

The required method:

* `Base.getindex(submat, x, y) = <substitution score/cost from x to y>`
"""
abstract AbstractSubstitutionMatrix{T<:Real}

function show_submat(io::IO, submat::AbstractSubstitutionMatrix, alphabet, sz)
    alphabets = [string(alphabet(UInt8(x))) for x in 0:sz-1]
    mat = [string(submat[x,y]) for x in 0:sz-1, y in 0:sz-1]
    # add rows
    mat = hcat(alphabets, mat)
    # add columns
    mat = vcat(vcat(" ", alphabets)', mat)
    width = [maximum([length(mat[i,j]) + 1 for j in 1:sz+1]) for i in 1:sz+1]
    println(io, typeof(submat), ":")
    for i in 1:sz+1
        print(io, "  ")  # indent
        for j in 1:sz+1
            el = mat[i,j]
            # right aligned
            print(io, " " ^ (width[j] - length(el)), el)
        end
        println(io)
    end
end


"""
Substitution matrix.
"""
type SubstitutionMatrix{T} <: AbstractSubstitutionMatrix{T}
    data::Matrix{T}
    alphabet::Nullable{DataType}
end

function Base.convert{T}(::Type{SubstitutionMatrix}, submat::AbstractMatrix{T})
    return SubstitutionMatrix{T}(submat, Nullable())
end

function SubstitutionMatrix{T}(submat::AbstractMatrix{T}, alphabet::DataType)
    return SubstitutionMatrix{T}(submat, Nullable(alphabet))
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

function Base.show(io::IO, submat::SubstitutionMatrix)
    alphabet = get(submat.alphabet, UInt8)
    sz = size(submat.data, 1)
    show_submat(io, submat, alphabet, sz)
end


"""
Dichotomous substitution matrix.
"""
immutable DichotomousSubstitutionMatrix{T} <: AbstractSubstitutionMatrix{T}
    match::T
    mismatch::T
end

function DichotomousSubstitutionMatrix(match::Real, mismatch::Real)
    typ = promote_type(typeof(match), typeof(mismatch))
    return DichotomousSubstitutionMatrix{typ}(match, mismatch)
end

@inline function Base.getindex(submat::DichotomousSubstitutionMatrix, x, y)
    return ifelse(
        x == y,
        submat.match,
        submat.mismatch
    )
end

function Base.show(io::IO, submat::DichotomousSubstitutionMatrix)
    match = string(submat.match)
    mismatch = string(submat.mismatch)
    w = max(length(match), length(mismatch))
    println(io, typeof(submat), ':')
    # right aligned
    print(io, "     match = ", " " ^ (w - length(match)), match, '\n')
    print(io, "  mismatch = ", " " ^ (w - length(mismatch)), mismatch)
end


# Utils for loading substitution matrices from files

function load_submat(name)
    submatfile = Pkg.dir("Bio", "src", "align", "data", "submat", name)
    return parse_ncbi_submat(submatfile)
end

function parse_ncbi_submat(filepath)
    ncols = 24
    # these amino acid code are not supported by Bio.Seq
    unsupported_aa = ('B', 'Z', '*')
    sz = ncols - length(unsupported_aa)
    submat = SubstitutionMatrix(Array{Int}(sz, sz), AminoAcid)
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
                if a in unsupported_aa
                    continue
                end
                for col in 1:ncols
                    b = header[col]
                    if b in unsupported_aa
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
