module Indexing

export
    Indexer, names!, rename!, rename


"""
The Indexer type

The type is parametric and contains a dictionary that maps
names (Symbols) to objects, as well as a vector of names
(Symbols) that enables indexers to have order (i.e.
The user should be able to push and pop and insert with indexes).
"""
immutable Indexer{T}
    lookup::Dict{Symbol, T}
    names::Vector{Symbol}
end


"Construct an empty Indexer that maps Symbols to a given Type T"
function Indexer{T}(::Type{T})
    return Indexer{T}(Dict{Symbol, T}(), Vector{Symbol}())
end


"""
Create an index that associates names in a Vector of Symbols, with values of type
T, in a Vector of type T.

This constructor constructs the index such that names[1] -> vals[1],
names[2] -> vals[2] ... and so on.
"""
function Indexer{T}(names::Vector{Symbol}, vals::Vector{T})
    u = make_unique(names)
    lookup = Dict{Symbol, T}(zip(u, vals))
    return Indexer{T}(lookup, u)
end


"""
Create an index that associates names in a Vector of Strings, with values of
type T, in a Vector of type T.

This constructor constructs the index such that names[1] -> vals[1],
names[2] -> vals[2] ... and so on.
"""
function Indexer{S <: AbstractString, T}(names::Vector{S}, vals::Vector{T})
    Indexer(convert(Vector{Symbol}, names), vals)
end


"""
Construct an Indexer that maps a Vector of Symbols (names), to a set of
values of type T <: Integer.

Provide this constructor a vector of Symbols as names, and a type that is a
subtype of Integer, and it will construct an index such that
names[1] -> 1, names[2] -> 2, ..., names[length(names)] -> length(names).
"""
function Indexer{T <: Integer}(names::Vector{Symbol}, ::Type{T})
    return Indexer(names, collect(T(1):T(length(names))))
end


"""
Construct an Indexer that maps a Vector of Strings (names), to a set of
values of type T <: Integer.

Provide this constructor a vector of Strings as names, and a type that is a
subtype of Integer, and it will construct an index such that
names[1] -> 1, names[2] -> 2, ..., names[length(names)] -> length(names).
"""
function Indexer{S <: AbstractString, T <: Integer}(names::Vector{S}, ::Type{T})
    return Indexer(convert(Vector{Symbol}, names), T)
end


"""
Construct an Indexer that maps a Vector of Symbols (names), to a set of
values of type Int, where Int is the default integer encoding on you machine
(usually Int32 or Int64).

Provide this constructor a vector of Symbols as names, and it will construct
an index such that names[1] -> 1, names[2] -> 2, ...,
names[length(names)] -> length(names).
"""
function Indexer(names::Vector{Symbol})
    return Indexer(names, Int)
end


"""
Construct an Indexer that maps a Vector of Strings (names), to a set of
values of type Int, where Int is the default integer encoding on you machine
(usually Int32 or Int64).

Provide this constructor a vector of Strings as names, and it will construct
an index such that names[1] -> 1, names[2] -> 2, ...,
names[length(names)] -> length(names).
"""
function Indexer{S <: AbstractString}(names::Vector{S})
    return Indexer(convert(Vector{Symbol}, names), Int)
end


# Internal function that makes sure all the symbols in a vector are unique.
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
    if length(x.names) != length(x.lookup)
        throw(ArgumentError("Index corrupted by duplicate symbols in names."))
    end
    rename!(x, x.names, names)
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
        nameidx = findfirst(x.names, from)
        x.names[nameidx] = to
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


# Indexing into an index mapping single names to single numbers

## Indexing with single name
Base.getindex(x::Indexer, idx::Symbol) = x.lookup[idx]
Base.getindex{S <: AbstractString}(x::Indexer, idx::S) = x.lookup[convert(Symbol, idx)]

## Indexing with multiple names.
Base.getindex(x::Indexer, idx::Vector{Symbol}) = [x.lookup[i] for i in idx]
function Base.getindex{S <: AbstractString}(x::Indexer, idx::Vector{S})
    return [x.lookup[convert(Symbol, i)] for i in idx]
end

## Indexing with a single integer
function Base.getindex{I <: Integer}(x::Indexer, idx::I)
    n = x.names[idx]
    return x[n]
end

## Indexing with a vector of intergers
function Base.getindex{I <: Integer}(x::Indexer, idx::Vector{I})
    ns = x.names[idx]
    return x[ns]
end

## Indexing with a range of integers
function Base.getindex{I <: Integer}(x::Indexer, idx::UnitRange{I})
    ns = x.names[idx]
    return x[ns]
end

## Indexing with a vector of bools.
function Base.getindex(x::Indexer, idx::Vector{Bool})
    ns = x.names[idx]
    return x[ns]
end

end
