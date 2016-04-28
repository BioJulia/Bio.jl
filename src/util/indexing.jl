module Indexing

export
    Indexer, names!, rename!, rename


intsAndVecs = Union{Unsigned,
                    AbstractVector{UInt8},
                    AbstractVector{UInt16},
                    AbstractVector{UInt32},
                    AbstractVector{UInt64},
                    AbstractVector{UInt128}}

# The indexer Type
type Indexer{T <: intsAndVecs}
    lookup::Dict{Symbol, T}
    names::Vector{Symbol}
end

# Conveinience constructors for creating Indexer that associate one name,
# with in integer value.
function Indexer{A <: AbstractVector, T <: Unsigned}(names::A{Symbol}, ::Type{T})
    u = make_unique(names)
    lookup = Dict{Symbol, Unsigned}(zip(u, UnitRange{T}(1, length(u))))
    Indexer{T}(lookup, u)
end
function Indexer{S <: AbstractString, T <: Unsigned}(names::Vector{S}, ::Type{T})
    Indexer(convert(Vector{Symbol}, names), T)
end

# Conveinience constructors for creating Indexer that associate one name,
# with several integer values.
function Indexer{A <: AbstractVector, E <: Unsigned}(names::Vector{Symbol}, groups::A{E})
    @assert length(names) == length(groups)
    u = make_unique(names)
    lookup = Dict{Symbol, A{E}}(zip(u, groups))
    Indexer(lookup, u)
end
function Indexer{S <: AbstractString, A <: AbstractVector, E <: Unsigned}(names::Vector{S}, groups::A{E})
    Indexer(convert(Vector{Symbol}, names), groups)
end
function Indexer{T <: intsAndVecs}(::Type{T})
    Indexer(Dict{Symbol, T}(), Symbol[])
end

function make_unique(names::Vector{Symbol})
    seen = Set{Symbol}()
    names = copy(names)
    dups = Int[]
    for i in 1:length(names)
        name = names[i]
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
    names = make_unique(names)
    if length(names) != length(x)
        throw(ArgumentError("Length of new names doesn't match length of Index."))
    end
    for i in 1:length(names)
        newname = names[i]
        oldname = x.names[i]
        if newname != oldname
            if x.lookup[oldname] >= i
                delete!(x.lookup, oldname)
            end
            x.lookup[newname] = i
            x.names[i] = newname
        end
    end
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
        x.lookup[to] = col = pop!(x.lookup, from)
        x.names[col] = to
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
