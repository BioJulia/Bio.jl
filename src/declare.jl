# Method Declarations
# ===================

# Method declaration macro.
macro declare(names)
    @assert names.head == :tuple
    base_names = Set(Base.names(Base))
    ex = Expr(:block)
    for name in names.args
        if name ∈ base_names
            error("method name '$(name)' is already used in Base")
        end
        push!(ex.args, :(function $(name) end))
        push!(ex.args, :(export $(name)))
    end
    return esc(ex)
end

# Base functions used in (and exported from) Bio.jl.
@declare (
    distance,
    seqname,
    hasseqname,
    sequence,
    hassequence,
    metadata,
    leftposition,
    hasleftposition,
    rightposition,
    hasrightposition,
    isoverlapping,
    metainfoval,
    metainfotag,
    header,
    isfilled
)
