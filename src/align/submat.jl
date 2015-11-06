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
    mat = [haskey(submat, x, y) ? string(submat[x,y]) : "NA" for x in 0:sz-1, y in 0:sz-1]
    # add rows
    mat = hcat(alphabets, mat)
    # add columns
    mat = vcat(vcat(" ", alphabets)', mat)
    width = [maximum([length(mat[i,j]) + 1 for j in 1:sz+1]) for i in 1:sz+1]
    println(io, typeof(submat), ":")
    for i in 1:sz+1
        print(io, " ")  # indent
        for j in 1:sz+1
            print(io, lpad(mat[i,j], width[j]))
        end
        println(io)
    end
end


"""
Substitution matrix.
"""
immutable SubstitutionMatrix{T} <: AbstractSubstitutionMatrix{T}
    # square substitution matrix
    data::Matrix{T}
    # a character is defined or not
    defined::BitVector
    # matching/mismatching score for undefined characters
    match::T
    mismatch::T
    # alphabet type (e.g. AminoAcid)
    alphabet::Nullable{DataType}

    function SubstitutionMatrix(data, defined, match, mismatch, alphabet=Nullable())
        sz = size(data, 1)
        @assert sz == size(data, 2)
        @assert sz == length(defined)

        # fill undefined cells with match/mismatch
        undefined = ~defined
        i = 0
        while (i = findnext(undefined, i + 1)) > 0
            for j in 1:sz
                score = i == j ? match : mismatch
                data[i,j] = data[j,i] = score
            end
        end

        return new(data, defined, match, mismatch, alphabet)
    end
end

function Base.convert{T}(::Type{SubstitutionMatrix}, submat::AbstractMatrix{T})
    sz = size(submat, 1)
    return SubstitutionMatrix{T}(submat, trues(sz), 0, 0, Nullable())
end

function SubstitutionMatrix{T}(submat::AbstractMatrix{T};
                               defined::BitVector=trues(size(submat, 1)),
                               match=T(0), mismatch=T(0), alphabet=Nullable())
    return SubstitutionMatrix{T}(submat, defined, match, mismatch, alphabet)
end

Base.convert(::Type{Matrix}, submat::SubstitutionMatrix) = submat.data
Base.convert{T}(::Type{Matrix{T}}, submat::SubstitutionMatrix{T}) = submat.data

@inline function Base.getindex(submat::SubstitutionMatrix, x, y)
    return submat.data[UInt8(x)+1,UInt8(y)+1]
end

function Base.haskey(submat::SubstitutionMatrix, x)
    return submat.defined[UInt8(x)+1]
end

function Base.haskey(submat::SubstitutionMatrix, x, y)
    return haskey(submat, x) && haskey(submat, y)
end

Base.minimum(submat::SubstitutionMatrix) = minimum(submat.data)
Base.maximum(submat::SubstitutionMatrix) = maximum(submat.data)

function Base.show(io::IO, submat::SubstitutionMatrix)
    alphabet = get(submat.alphabet, UInt8)
    show_submat(io, submat, alphabet, size(submat.data, 1))
    println(io, "* NA: match = ", submat.match, ", mismatch = ", submat.mismatch)
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
    d = Dict{Tuple{Char,Char},Int}()
    # column => one-letter amino acid
    header = Char[]
    open(filepath) do io
        # column width (inferred from the header line)
        w = 0
        local ncols::Int
        for line in eachline(io)
            if startswith(line, "#")
                continue
            elseif isempty(header)
                # header line
                ncols = length(matchall(r"[A-Z*]", line))
                w = div(length(line) - 2, ncols)
                for col in 1:ncols
                    char = line[w*col+1]
                    @assert isalpha(char) || char == '*'
                    push!(header, char)
                end
            else
                a = line[1]
                for col in 1:ncols
                    b = header[col]
                    d[(a, b)] = parse(Int, line[w*col-1:w*col+1])
                end
            end
        end
    end
    # create the substitution matrix
    aas = alphabet(AminoAcid)
    n_aas = length(aas)
    defined = falses(n_aas)
    for char in header
        if char != '*'
            defined[UInt8(AminoAcid(char))+1] = true
        end
    end
    submat = Matrix{Int}(n_aas, n_aas)
    for (i, x) in enumerate(aas), (j, y) in enumerate(aas)
        submat[i,j] = get(d, (Char(x), Char(y)), 0)
    end
    return SubstitutionMatrix(
        submat,
        defined=defined,
        match=0,
        mismatch=0,
        alphabet=AminoAcid
    )
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
