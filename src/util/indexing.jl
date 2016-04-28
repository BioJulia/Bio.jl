module Indexing

export
    AbstractIndex, SingleIndex, names!, rename!, rename


# Indexing datastructures

# Abstract Types
abstract AbstractIndex


# Concrete Types
type SingleIndex <: AbstractIndex
    lookup::Dict{Symbol, Int}
    names::Vector{Symbol}
end

function SingleIndex(names::Vector{Symbol})
    u = make_unique(names)
    lookup = Dict{Symbol, Int}(zip(u, 1:length(u)))
    SingleIndex(lookup, u)
end
SingleIndex{S <: AbstractString}(names::Vector{S}) = SingleIndex(convert(Vector{Symbol}, names))
SingleIndex() = SingleIndex(Dict{Symbol, Int}(), Symbol[])


type GroupIndex{T <: AbstractVector{Int}} <: AbstractIndex
    lookup::Dict{Symbol, Int}
    names::Vector{Symbol}
    groups::Vector{T}
end

function GroupIndex{T <: AbstractVector{Int}}(names::Vector{Symbol}, groups::Vector{T})
    u = make_unique(names)
    lookup = Dict{Symbol, Int}(zip(u, 1:length(u)))
    GroupIndex(lookup, u, groups)
end
GroupIndex{S <: AbstractString, T <: AbstractVector{Int}}(names::Vector{S}, groups::Vector{T}) = GroupIndex(convert(Vector{Symbol}, names), groups)
GroupIndex{T <: AbstractVector{Int}}(::Type{T}) = GroupIndex(Dict{Symbol, Int}(), Symbol[], T[])

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

Base.length(x::AbstractIndex) = length(x.names)
Base.names(x::AbstractIndex) = copy(x.names)
_names(x::AbstractIndex) = x.names
Base.copy(x::SingleIndex) = SingleIndex(copy(x.lookup), copy(x.names))
Base.copy(x::GroupIndex) = GroupIndex(copy(x.lookup), copy(x.names), copy(x.groups))
Base.deepcopy(x::AbstractIndex) = copy(x) # all eltypes immutable
Base.isequal(x::SingleIndex, y::SingleIndex) = isequal(x.lookup, y.lookup) && isequal(x.names, y.names)
Base.isequal(x::GroupIndex, y::GroupIndex) = isequal(x.lookup, y.lookup) && isequal(x.names, y.names) && isequal(x.groups, y.groups)
Base.isequal(x::SingleIndex, y::GroupIndex) = false
Base.isequal(x::GroupIndex, y::SingleIndex) = false
Base.(:(==))(x::AbstractIndex, y::AbstractIndex) = isequal(x, y)
Base.haskey(x::AbstractIndex, key::Symbol) = haskey(x.lookup, key)
Base.haskey{T <: AbstractString}(x::AbstractIndex, key::T) = haskey(x, convert(Symbol, key))
Base.haskey(x::AbstractIndex, key::Real) = 1 <= key <= length(x.names)
Base.keys(x::AbstractIndex) = names(x)

function names!(x::AbstractIndex, names::Vector{Symbol})
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

function names!{T <: AbstractString}(x::AbstractIndex, names::Vector{T})
    names!(x, convert(Vector{Symbol}, names))
end

function rename!(x::AbstractIndex, names)
    for (from, to) in names
        if haskey(x, to)
            error("Tried renaming $from to $to, when $to already exists in the Index.")
        end
        x.lookup[to] = col = pop!(x.lookup, from)
        x.names[col] = to
    end
    return x
end

function rename!(x::AbstractIndex, from::Vector{Symbol}, to::Vector{Symbol})
    return rename!(x, zip(from, to))
end

rename!(x::AbstractIndex, from::Symbol, to::Symbol) = rename!(x, ((from, to),))
rename!(x::AbstractIndex, f::Function) = rename!(x, [(x,f(x)) for x in x.names])
rename!(f::Function, x::AbstractIndex) = rename!(x, f)
rename(x::AbstractIndex, args...) = rename!(copy(x), args...)
rename(f::Function, x::AbstractIndex) = rename(copy(x), f)


# Indexing into a single index.

# Indexing with single name.
Base.getindex(x::SingleIndex, idx::Symbol) = x.lookup[idx]
Base.getindex{S <: AbstractString}(x::SingleIndex, idx::S) = x.lookup[convert(Symbol, idx)]
# Indexing with a single integer.
Base.getindex(x::SingleIndex, idx::Int) = idx
# Indexing with a vector of bools.
Base.getindex(x::SingleIndex, idx::Vector{Bool}) = find(idx)
# Indexing with a range.
Base.getindex(x::SingleIndex, idx::UnitRange{Int}) = [idx;]
# Indexing with multiple integers.
Base.getindex{T <: Real}(x::SingleIndex, idx::Vector{T}) = convert(Vector{Int}, idx)
# Indexing with multiple names.
Base.getindex(x::SingleIndex, idx::Vector{Symbol}) = [x.lookup[i] for i in idx]
function Base.getindex{T <: AbstractString}(x::SingleIndex, idx::Vector{T})
    return [x.lookup[convert(Symbol, i)] for i in idx]
end


# Indexing into a grouped index.

# Indexing with a single name.
Base.getindex(x::GroupIndex, idx::Symbol) = x.groups[x.lookup[idx]]
Base.getindex{S <: AbstractString}(x::GroupIndex, idx::S) = x.groups[x.lookup[convert(Symbol, idx)]]
# Indexing with a single integer.
Base.getindex(x::GroupIndex, idx::Int) = x.groups[idx]
# Indexing with a vector of bools.
Base.getindex(x::GroupIndex, idx::Vector{Bool}) = x.groups[find(idx)]
# Indexing with a range.
Base.getindex(x::GroupIndex, idx::Range) = x.groups[[idx;]]
# Indexing with multiple integers.
Base.getindex{T <: Real}(x::GroupIndex, idx::Vector{T}) = x.groups[convert(Vector{Int}, idx)]
# Indexing with multiple names.
Base.getindex(x::GroupIndex, idx::Vector{Symbol}) = [x.groups[x.lookup[i]] for i in idx]

end
