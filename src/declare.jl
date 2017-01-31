# Method Declarations
# ===================

# Method declaration macro.
macro declare(names)
    @assert names.head == :tuple
    ex = Expr(:block)
    for name in names.args
        push!(ex.args, :(function $(name) end))
        push!(ex.args, :(export $(name)))
    end
    return esc(ex)
end

@declare (
    distance,
)
