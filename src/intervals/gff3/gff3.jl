
# Generic Feature Format Version 3 (GFF3)
# =======================================
#
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

using URIParser: unescape

"""
Map strings to one or more strings avoiding allocation and deallocation as much
as possible.
"""
type GFFAttributes <:
        Associative{StringField, Union{StringField, Vector{StringField}}}
    # map key to corresponding index in data and used
    indexes::Dict{StringField, Int}

    data::Vector{Vector{StringField}}

    # how many fields in data[i] are used
    used::Vector{Int}

    function GFFAttributes()
        return new(Dict{StringField, Int}(), Vector{StringField}[], Int[])
    end
end

function Base.length(attrs::GFFAttributes)
    return sum(attrs.used)
end

function Base.empty!(attrs::GFFAttributes)
    fill!(attrs.used, 0)
end

function Base.haskey(attrs::GFFAttributes, key::StringField)
    i = get(attrs.indexes, key, 0)
    return i != 0 && attrs.used[i] > 0
end

function Base.haskey(attrs::GFFAttributes, key_::String)
    key = StringField(key_.data, 1:length(key_.data))
    return haskey(attrs, key)
end

function Base.getindex(attrs::GFFAttributes, key::StringField)
    i = get(attrs.indexes, key, 0)
    if i == 0 || !attrs.used[i]
        throw(KeyError(key))
    end
    return attrs.data[i]
end

function Base.getindex(attrs::GFFAttributes, key_::String)
    key = StringField(key_.data, 1:length(key_.data))
    return getindex(attrs, key)
end

function pushindex!(attrs::GFFAttributes, key::StringField,
                    data::Vector{UInt8}, start::Int, stop::Int,
                    unescape_needed::Bool)
    i = get(attrs.indexes, key, 0)
    if i == 0
        i = attrs.indexes[key] = length(attrs.indexes) + 1
        push!(attrs.used, 0)
        push!(attrs.data, StringField[])
    end
    j = (attrs.used[i] += 1)
    if j > length(attrs.data[i])
        push!(attrs.data[i], StringField())
    end
    copy!(attrs.data[i][j], data, start, stop)
    if unescape_needed
        unescape!(attrs.data[i][j])
    end
end

function Base.copy(attrs::GFFAttributes)
    attrs2 = GFFAttributes()
    for (key, value) in attrs
        i = get(attrs2.indexes, key, 0)
        if i == 0
            i = attrs2.indexes[key] = length(attrs2.indexes) + 1
            push!(attrs2.used, 0)
            push!(attrs2.data, StringField[])
        end
        j = (attrs2.used[i] += 1)
        push!(attrs2.data[i], copy(value))
    end
    return attrs2
end

function Base.start(attrs::GFFAttributes)
    index_iter_state = start(attrs.indexes)
    i, j, key = 1, 1, StringField()
    while !done(attrs.indexes, index_iter_state) && (i == 0 || j > attrs.used[i])
        j = 1
        keyindex, index_iter_state = next(attrs.indexes, index_iter_state)
        key, i = keyindex.first, keyindex.second
    end
    return (index_iter_state, key, i, j)
end

function Base.next(attrs::GFFAttributes, state)
    index_iter_state, key, i, j = state
    value = attrs.data[i][j]
    j += 1
    while j > attrs.used[i] && !done(attrs.indexes, index_iter_state)
        j = 1
        keyindex, index_iter_state = next(attrs.indexes, index_iter_state)
        key, i = keyindex.first, keyindex.second
    end
    return Pair{StringField, StringField}(key, value), (index_iter_state, key, i, j)
end

function Base.done(attrs::GFFAttributes, state)
    index_iter_state, key, i, j = state
    return i > length(attrs.used) ||
        (j > attrs.used[i] && done(attrs.indexes, index_iter_state))
end

type GFF3Metadata
    source::StringField
    # use "kind" instead of "type" since "type" is a keyword of Julia
    kind::StringField
    score::Nullable{Float64}
    phase::Nullable{Int}
    attributes::GFFAttributes
end

function GFF3Metadata()
    return GFF3Metadata("", "", NaN, 0, GFFAttributes())
end

function Base.copy(metadata::GFF3Metadata)
    return GFF3Metadata(
        copy(metadata.source),
        copy(metadata.kind),
        metadata.score,
        metadata.phase,
        copy(metadata.attributes))
end

function unescape!(x::StringField)
    y = unescape(x)
    copy!(x, y.data, 1, length(y.data))
    return x
end

"An `Interval` with associated metadata from a GFF3 file"
typealias GFF3Interval Interval{GFF3Metadata}

include("reader.jl")
include("parser.jl")
