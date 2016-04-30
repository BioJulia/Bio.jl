module Indexing

export
    Indexer, names!, rename!, rename


# A bit of a clunky workaround until we get triangle dispatch in 0.5 or 0.6
indexerUnion = Union{Unsigned, AbstractVector{UInt8}, AbstractVector{UInt16},
                     AbstractVector{UInt32}, AbstractVector{UInt64},
                     AbstractVector{UInt128}}

vectorsUnion = Union{AbstractVector{UInt8}, AbstractVector{UInt16},
                     AbstractVector{UInt32}, AbstractVector{UInt64},
                     AbstractVector{UInt128}}


# The indexer Type
immutable Indexer{T <: indexerUnion}
    lookup::Dict{Symbol, T}
    names::Vector{Symbol}
end

# The reason the indexer has a names field as well as a Dict is because
# we can concieve of cases in which the indexer could be considered ordered,
# and so operations like pushing, poping, and inserting may be applicable.
# mostly in cases where the indexer is used with genotype matrices.

function Indexer{T <: indexerUnion}(::Type{T})
    return Indexer{T}(Dict{Symbol, T}(), Vector{Symbol}())
end

# Convenience constructors for creating Indexer that associate one name,
# with in integer value.
function Indexer{T <: Unsigned}(names::Vector{Symbol}, ::Type{T})
    u = make_unique(names)
    lookup = Dict{Symbol, T}(zip(u, UnitRange{T}(1, length(u))))
    return Indexer{T}(lookup, u)
end

function Indexer(names::Vector{Symbol})
    return Indexer(names, UInt)
end

function Indexer{S <: AbstractString}(names::Vector{S})
    return Indexer(convert(Vector{Symbol}, names), UInt)
end

function Indexer{S <: AbstractString, T <: Unsigned}(names::Vector{S}, ::Type{T})
    return Indexer(convert(Vector{Symbol}, names), T)
end

# Conveinience constructors for creating Indexer that associate one name,
# with several integer values.
function Indexer{T <: vectorsUnion}(names::Vector{Symbol}, groups::Vector{T})
    @assert length(names) == length(groups)
    u = make_unique(names)
    lookup = Dict{Symbol, T}(zip(u, groups))
    Indexer{T}(lookup, u)
end

function Indexer{S <: AbstractString, T <: vectorsUnion}(names::Vector{S}, groups::Vector{T})
    Indexer(convert(Vector{Symbol}, names), groups)
end

function make_unique(names::Vector{Symbol})
    seen = Set{Symbol}()
    names = copy(names)
    dups = Int[]
    for (i, name) in enumerate(names)
        in(name, seen) ? push!(dups, i) : push!(seen, name)
    end
    for i in dups
        nm = names[i]
        k = 1
        while true
            newnm = symbol("$(nm)_$k")
            if !in(newnm, seen)
                warn("Trying to build an index with duplicate names. Changing $(names[i]) to $(newnm)")
                names[i] = newnm
                push!(seen, newnm)
                break
            end
            k += 1
        end
    end
    return names
end

Base.length(x::Indexer) = length(x.names)
Base.names(x::Indexer) = copy(x.names)
_names(x::Indexer) = x.names
Base.copy(x::Indexer) = Indexer(copy(x.lookup), copy(x.names))

Base.deepcopy(x::Indexer) = copy(x) # all eltypes immutable
Base.isequal(x::Indexer, y::Indexer) = isequal(x.lookup, y.lookup) && isequal(x.names, y.names)
Base.(:(==))(x::Indexer, y::Indexer) = isequal(x, y)
Base.haskey(x::Indexer, key::Symbol) = haskey(x.lookup, key)
Base.haskey{S <: AbstractString}(x::Indexer, key::S) = haskey(x, convert(Symbol, key))
Base.haskey(x::Indexer, key::Real) = 1 <= key <= length(x.names)
Base.keys(x::Indexer) = names(x)

function names!(x::Indexer, names::Vector{Symbol})
    if length(names) != length(x)
        throw(ArgumentError("Length of new names doesn't match length of Index."))
    end
    names = make_unique(names)
    rename!(x, names, x.names)
    if length(x.names) != length(x.lookup)
        throw(ArgumentError("Index corrupted by duplicate symbols in names."))
    end
    return x
end

function names!{S <: AbstractString}(x::Indexer, names::Vector{S})
    names!(x, convert(Vector{Symbol}, names))
end

function rename!(x::Indexer, names)
    for (from, to) in names
        if haskey(x, to)
            error("Tried renaming $from to $to, when $to already exists in the Index.")
        end
        # Change the name array
        nameidx = findin(x.names, from)
        names[nameidx] = to
        # Change the dictionary
        x.lookup[to] = val = pop!(x.lookup, from)
    end
    return x
end

function rename!(x::Indexer, from::Vector{Symbol}, to::Vector{Symbol})
    return rename!(x, zip(from, to))
end

rename!(x::Indexer, from::Symbol, to::Symbol) = rename!(x, ((from, to),))
rename!(x::Indexer, f::Function) = rename!(x, [(x,f(x)) for x in x.names])
rename!(f::Function, x::Indexer) = rename!(x, f)
rename(x::Indexer, args...) = rename!(copy(x), args...)
rename(f::Function, x::Indexer) = rename(copy(x), f)


# Indexing into a single index.

# Indexing with single name.
Base.getindex(x::Indexer, idx::Symbol) = x.lookup[idx]
Base.getindex{S <: AbstractString}(x::Indexer, idx::S) = x.lookup[convert(Symbol, idx)]
# Indexing with a single integer.
Base.getindex{T <: Unsigned}(x::Indexer{T}, idx::T) = idx
# Indexing with a vector of bools.
Base.getindex{T <: Unsigned}(x::Indexer{T}, idx::Vector{Bool}) = T(find(idx))
# Indexing with a range.
Base.getindex{T <: Unsigned}(x::Indexer, idx::UnitRange{T}) = [idx;]
# Indexing with multiple names.
Base.getindex(x::Indexer, idx::Vector{Symbol}) = [x.lookup[i] for i in idx]
function Base.getindex{S <: AbstractString}(x::Indexer, idx::Vector{S})
    return [x.lookup[convert(Symbol, i)] for i in idx]
end


end
