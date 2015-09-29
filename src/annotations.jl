# Annotations
# ===========

#=
Type-stable, yet flexible container. Useful for storing things like phylogeny
and sequence annotations, without hving to simply accept the speed hit or type
instability.

Design borrows heavily from DataFrames prototypes for type stable DataFrames.
Because the type of a tuple is a tupleof the types, we can exploit this to get
some type stability and predictability for a flexible data structure in which
types may not be known in advance.
=#

parameters{T}(::Type{T}) = T.parameters

abstract AbstractAnnotation

type NoAnnotations <: AbstractAnnotation
end

type Annotations{D <: Tuple, L <: Tuple} <: AbstractAnnotation
    data::D
    labels::Type{L}
end


# Generated function outer constructor prevents duplicate labels for fields.
@generated function Annotations{D <: Tuple, L <: Tuple}(data::D, labels::Type{L})
    lparams = parameters(L)
    if length(unique(lparams)) < length(lparams)
        return :(error("Trying to give more than one field the same label."))
    else
        return :(Annotations{D, L}(data, labels))
    end
end

# Type stable construction is made much more easy with a macro:

macro annotations(kwargs...)
    labels = map(i->i.args[1], kwargs)
    data = Expr(:tuple, map(i->i.args[2], kwargs)...)
    return :(Annotations($data, Tuple{$(labels...)}))
end

abstract Field{a}

@generated function Base.getindex{D, L, a}(annotation::Annotations{D, L}, ::Type{Field{a}})
    index = findfirst(parameters(L), a)
    return :(annotation.data[$index])
end

Base.getindex{a}(annotations::Annotations, i, ann::Type{a}) = annotations[field][i]

@generated function Base.setindex!{D, L, A, a}(annotations::Annotations{D, L}, addition::A, ::Type{Field{a}})
    L2 = Tuple{parameters(L)..., a}
    D2 = Tuple{parameters(D)..., A}
    return :(Annotations{$D2, $L2}(tuple(annotations.data..., addition), $L2))
end
