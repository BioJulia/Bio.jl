# Indexers
# ========
#
# Mapping from names to objects.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Indexers

export
    Indexer, names!, rename!, rename

using Compat

"""
The Indexer type

The type is parametric and contains a dictionary mapping
names (Symbols) to objects, as well as a vector of names
(Symbols) that enables indexers to have order (i.e.
The user should be able to push and pop and insert with indexes).
"""
immutable Indexer{T}
    lookup::Dict{Symbol, T}
    names::Vector{Symbol}
end


"""
    Indexer{T}(::Type{T})

Construct an empty Indexer mapping Symbols to a given Type T.

# Examples
```julia
Indexer(Int64)
Indexer(Int8)
Indexer(UInt32)
```
"""
function Indexer{T}(::Type{T})
    return Indexer{T}(Dict{Symbol, T}(), Vector{Symbol}())
end


"""
    Indexer{T}(names::Vector{Symbol}, vals::Vector{T})

Create an index that associates a Vector of Symbols with a vector of values.

This constructor constructs the index such that names[1] -> vals[1],
names[2] -> vals[2] ... and so on.

# Examples
```julia
Indexer([:Hello, :Bye], [1, 2])
Indexer([:Hello, :Bye], [1:10, 11:20])
```
"""
function Indexer{T}(names::Vector{Symbol}, vals::Vector{T})
    u = make_unique(names)
    lookup = Dict{Symbol, T}(zip(u, vals))
    return Indexer{T}(lookup, u)
end


"""
    Indexer{S <: AbstractString, T}(names::Vector{S}, vals::Vector{T})

Create an index that associates names in a Vector of Strings with a vector of
values.

This constructor constructs the index such that names[1] -> vals[1],
names[2] -> vals[2] ... and so on.

# Examples
```julia
Indexer(["Hi", "Bye", "Sup"], [1, 2, 3])
Indexer(["Hi", "Bye", "Sup"], [1:5, 6:10, 11:15])
```
"""
function Indexer{S <: AbstractString, T}(names::Vector{S}, vals::Vector{T})
    Indexer(convert(Vector{Symbol}, names), vals)
end


"""
    Indexer(names::Vector{Symbol}, ::Type{T})

Construct an Indexer mapping a Vector of Symbols (names), to a set of
values of type T.

Provide this constructor a vector of Symbols as names, and a type that is a
subtype of Integer, and it will construct an index such that
names[1] -> 1, names[2] -> 2, ..., names[length(names)] -> length(names).

# Examples
```julia
Indexer([:First, :Second, :Third], Int64)
Indexer([:First, :Second, :Third], UInt32)
```
"""
function Indexer{T <: Integer}(names::Vector{Symbol}, ::Type{T})
    return Indexer(names, collect(T(1):T(length(names))))
end


"""
    Indexer{S <: AbstractString, T <: Integer}(names::Vector{S}, ::Type{T})

Construct an Indexer mapping a Vector of Strings (names), to a set of
values of type T <: Integer.

Provide this constructor a vector of Strings as names, and a type that is a
subtype of Integer, and it will construct an index such that
names[1] -> 1, names[2] -> 2, ..., names[length(names)] -> length(names).

# Examples
```julia
Indexer(["First", "Second", "Third"], Int64)
Indexer(["First", "Second", "Third"], UInt32)
```
"""
function Indexer{S <: AbstractString, T <: Integer}(names::Vector{S}, ::Type{T})
    return Indexer(convert(Vector{Symbol}, names), T)
end


"""
    Indexer(names::Vector{Symbol})

Construct an Indexer mapping a Vector of Symbols (names), to a set of
values of type Int, where Int is the default integer encoding on your machine
(usually Int32 or Int64).

Provide this constructor a vector of Symbols as names, and it will construct
an index such that names[1] -> 1, names[2] -> 2, ...,
names[length(names)] -> length(names).

# Examples
```julia
Indexer([:First, :Second, :Third])
```
"""
function Indexer(names::Vector{Symbol})
    return Indexer(names, Int)
end


"""
    Indexer{S <: AbstractString}(names::Vector{S})

Construct an Indexer mapping a Vector of Strings (names), to a set of
values of type Int, where Int is the default integer encoding on you machine
(usually Int32 or Int64).

Provide this constructor a vector of Strings as names, and it will construct
an index such that names[1] -> 1, names[2] -> 2, ...,
names[length(names)] -> length(names).

# Examples
```julia
Indexer(["First", "Second", "Third"])
```
"""
function Indexer{S <: AbstractString}(names::Vector{S})
    return Indexer(convert(Vector{Symbol}, names), Int)
end


"""
    make_unique(names::Vector{Symbol})

Make sure all vectors in a symbol are unique.

# Examples
```julia
make_unique([:First, :Second, :Third])
make_unique([:First, :First, :First])
make_unique([:First, :Third, :Third])
```
"""
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
            newnm = Symbol("$(nm)_$k")
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
Base.copy(x::Indexer) = Indexer(copy(x.lookup), copy(x.names))
@compat Base.:(==)(x::Indexer, y::Indexer) = (x.lookup == y.lookup) && (x.names == y.names)
Base.haskey(x::Indexer, key::Symbol) = haskey(x.lookup, key)
Base.haskey(x::Indexer, key::AbstractString) = haskey(x, convert(Symbol, key))
Base.haskey(x::Indexer, key::Integer) = 1 <= key <= length(x.names)
Base.keys(x::Indexer) = names(x)


"""
    names!(x::Indexer, names::Vector{Symbol})

Completely replace the names in an Indexer with a new Vector of Symbols.

# Examples
```julia
names!(my_indexer, [:First, :NewName, :Third])
```
"""
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


"""
    names!{S <: AbstractString}(x::Indexer, names::Vector{S})

Completely replace the names in an Indexer with a new Vector of Strings.

# Examples
```julia
names!(my_indexer, ["First", "NewName", "Third"])
```
"""
function names!{S <: AbstractString}(x::Indexer, names::Vector{S})
    names!(x, convert(Vector{Symbol}, names))
end


function rename!(x::Indexer, names)
    for (from, to) in names
        if !haskey(x, from)
            throw(ArgumentError("Cannot rename $from to $to, $from does not exist in the Indexer."))
        end
        if haskey(x, to)
            throw(ArgumentError("Cannot rename $from to $to, $to already exists in the Indexer."))
        end
        # Change the name array
        nameidx = findfirst(x.names, from)
        x.names[nameidx] = to
        # Change the dictionary
        x.lookup[to] = val = pop!(x.lookup, from)
    end
    return x
end


"""
    rename!(x::Indexer, from::Vector{Symbol}, to::Vector{Symbol})

Rename a set of names in an Indexer, to a new set of names.

This will rename such that the name in from[1] will be renamed to to[1],
from[2] will be renamed to to[2] and so on.

# Examples
```julia
rename!(my_indexer, [:First, :Third], [:NewFirst, :NewThird])
```
"""
function rename!(x::Indexer, from::Vector{Symbol}, to::Vector{Symbol})
    return rename!(x, zip(from, to))
end


"""
    rename!{S <: AbstractString}(x::Indexer, from::Vector{S}, to::Vector{S})

Rename a set of names in an Indexer, to a new set of names.

This will rename such that the name in from[1] will be renamed to to[1],
from[2] will be renamed to to[2] and so on.

# Examples
```julia
rename!(my_indexer, ["First", "Third"], ["NewFirst", "NewThird"])
```
"""
function rename!{S <: AbstractString}(x::Indexer, from::Vector{S}, to::Vector{S})
    return rename!(x, convert(Vector{Symbol}, from), convert(Vector{Symbol}, to))
end


"""
    rename!(x::Indexer, from::Symbol, to::Symbol)

Rename a single name (from) in Indexer x, to a new name (to).

# Examples
```julia
rename!(my_indexer, :First, :NewFirst)
```
"""
rename!(x::Indexer, from::Symbol, to::Symbol) = rename!(x, ((from, to),))

"""
    rename!{S <: AbstractString}(x::Indexer, from::S, to::S)

Rename a single name (from) in Indexer x, to a new name (to).

# Examples
```julia
rename!(my_indexer, "First", "NewFirst")
```
"""
function rename!{S <: AbstractString}(x::Indexer, from::S, to::S)
    return rename!(x, ((convert(Symbol, from), convert(Symbol, to)),))
end


"""
    rename!(f::Function, x::Indexer)

Rename a indexer according to some function, f, that accepts one argument.

Each name in the Indexer x, will be input to, f, and the output will be
what that name will be renamed to.
Therefore the function, f, will need to return Symbols or Strings.
"""
rename!(f::Function, x::Indexer) = rename!(x, [(x,f(x)) for x in x.names])
rename(x::Indexer, args...) = rename!(copy(x), args...)


# Indexing into an index mapping single names to single numbers

## Indexing with single name
Base.getindex(x::Indexer, key::Symbol) = x.lookup[key]
Base.getindex(x::Indexer, key::AbstractString) = x.lookup[convert(Symbol, key)]

## Indexing with multiple names.
Base.getindex(x::Indexer, keys::Vector{Symbol}) = [x.lookup[key] for key in keys]
function Base.getindex{S <: AbstractString}(x::Indexer, keys::Vector{S})
    return [x.lookup[convert(Symbol, key)] for key in keys]
end

## Indexing with a single integer
function Base.getindex{I <: Integer}(x::Indexer, key::I)
    n = x.names[key]
    return x[n]
end

## Indexing with a range of integers
function Base.getindex{I <: Integer}(x::Indexer, keys::AbstractVector{I})
    ns = x.names[keys]
    return x[ns]
end

## Indexing with a vector of bools.
function Base.getindex(x::Indexer, keys::Vector{Bool})
    ns = x.names[keys]
    return x[ns]
end

end
