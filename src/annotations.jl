# Design borrows heavily from DataFrames prototypes for type stable DataFrames.
# Because the type of a tuple is a tupleof the types, we can exploit this to get
# some type stability and predictability for a flexible data structure in which types
# may not be known in advance.

parameters{T}(::Type{T}) = T.parameters

abstract AbstractAnnotation

type NoAnnotations <: AbstractAnnotation
end

type Annotations{D <: Tuple, L <: Tuple} <: AbstractAnnotation
    data::D
    labels::Type{L}
end

# Type stable construction is possible with a macro:

macro annotations(kwargs...)
    labels = map(i->i.args[1], kwargs)
    data = Expr(:tuple, map(i->i.args[2], kwargs)...)
    return :(Annotations($data, Tuple{$(labels...)}))
end

abstract Annotation{a}

@generated function Base.getindex{D, L, a}(annotation::Annotations{D, L}, ::Type{Annotation{a}})
    index = findfirst(parameters(L), a)
    return :(annotation.data[$index])
end

Base.getindex{a}(annotations::Annotations, i, ann::Type{a}) = annotations[field][i]

@generated function Base.setindex!{D, L, A, a}(annotations::Annotations{D, L}, addition::A, ::Type{Annotation{a}})
    L2 = Tuple{parameters(L)..., a}
    D2 = Tuple{parameters(D)..., A}
    return :(Annotations{$D2, $L2}(tuple(annotations.data..., addition), $L2))
end
