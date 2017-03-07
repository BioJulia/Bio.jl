# Exceptions
# ==========

module Exceptions

"""
    MissingFieldException <: Exception

An exception type thrown when a missing field of a record is accessed.
"""
immutable MissingFieldException <: Exception
    field::Symbol
end

# This shows a warning message?
function Base.showerror(io::IO, ex::MissingFieldException)
    print(io, "$(ex.field) is missing")
end

# Throws a MissingFieldException exception.
function missingerror(field::Symbol)
    throw(MissingFieldException(field))
end

end
