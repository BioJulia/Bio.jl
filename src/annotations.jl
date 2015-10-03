# Annotations
# ===========

abstract AbstractAnnotation

"""
Type-stable, yet flexible container. Useful for storing things like phylogeny
and sequence annotations, without having to simply accept the speed hit or type
instability that may occur.

Design borrows heavily from DataFrames prototypes for type-stable DataFrames.
Because the type of a tuple is a tuple of the types, we can exploit this to get
some type stability and predictability for a flexible data structure in which
types may not be known in advance.
"""
immutable Annotations{D <: Tuple, L <: Tuple} <: AbstractAnnotation
    "Stores a tuple of variables of any type.
    The type of this field then is a tuple of the types."
    data::D
    "Stores the labels used to getindex the data stored in the data variable"
    labels::Type{L}
end

"A convienient empty type, which inherits from the same abstract type as Annotations.
Used in other types to denote that no annotations are stored."
type NoAnnotations <: AbstractAnnotation end

"Abstract type used to subset into annotations."
abstract Field{a}

"Get the parameters of a given type."
parameters{T}(::Type{T}) = T.parameters

# Generated function outer constructor prevents duplicate labels for fields.
# Not sure if this needs to be a generated function by nessecity, but I don't
# think it hurts to move the labels checking to compile-time, given the type
# info is available at that point.
"""
Annotations(data, labels::Type{L})

Construct an Annotations container.

**Parameters:**

*   *data*

    The actual variables to be stored in the annotations object.
    They must be provided as a tuple.

*   *labels*

    Labels with which to identify the variables stored in the annotations object.
    This needs to be given as a Tuple Type, usually with symbols as parameters.
    For example: Tuple{:a, :b} for an Annotations object containing two variables.
    The two variables are assigned to fields :a and :b, or rather Field{:a} and Field{:b}.
"""
@generated function Annotations{D <: Tuple, L <: Tuple}(data::D, labels::Type{L})
    lparams = parameters(L)
    if length(unique(lparams)) < length(lparams)
        return :(error("Trying to give more than one field the same label."))
    else
        return :(Annotations{D, L}(data, labels))
    end
end

"Calling Annotations() without any arguments generates a NoAnnotations object."
Annotations() = NoAnnotations()

"Generated function outer constructor for extending annotations in a functional programming manner."
@generated function Annotations{D, L, A, a}(annotations::Annotations{D, L}, addition::A, ::Type{Field{a}})
    L2 = Tuple{parameters(L)..., a}
    D2 = Tuple{parameters(D)..., A}
    return :(Annotations{$D2, $L2}(tuple(annotations.data..., addition), $L2))
end



# Type stable construction is made much more easy with a macro:

"annotations Macro, allows easy construction of Annotations objects by providing a list of keyword-value pairs."
macro annotations(kwargs...)
    labels = map(i -> i.args[1], kwargs)
    data = Expr(:tuple, map(i -> i.args[2], kwargs)...)
    return :(Annotations($data, Tuple{$(labels...)}))
end


"Get any variable from an annotation variable field, whilst achieving type stability."
@generated function Base.getindex{D, L, a}(annotation::Annotations{D, L}, ::Type{Field{a}})
    index = findfirst(parameters(L), a)
    return :(annotation.data[$index])
end

"Get any variable from an annotation variable field, whilst achieving type stability."
Base.getindex{a}(annotations::Annotations, i, ann::Type{a}) = annotations[field][i]
