# SAM Header
# ==========

# Implicit type restriction for predefiend tags:
#   * HD: Dict{String,String}
#   * SQ: Vector{Dict{String,Any}}
#   * RG: Vector{Dict{String,Any}}
#   * PG: Vector{Dict{String,Any}}
#   * CO: Vector{String}

type SAMHeader <: Associative{String,Any}
    data::Dict{String,Any}
end

function SAMHeader()
    return SAMHeader(Dict{String,Any}())
end

function Base.length(header::SAMHeader)
    return length(header.data)
end

function Base.keys(header::SAMHeader)
    return keys(header.data)
end

function Base.values(header::SAMHeader)
    return values(header.data)
end

function Base.getindex(header::SAMHeader, key::AbstractString)
    return getindex(header.data, key)
end

function Base.setindex!(header::SAMHeader, val, key::AbstractString)
    if length(key) != 2
        error("SAM header tag must be of length 2")
    end
    # TODO: check value type and validity?
    setindex!(header.data, val, key)
    return header
end

function Base.delete!(header::SAMHeader, key::AbstractString)
    return delete!(header.data, key)
end

function Base.start(header::SAMHeader)
    return start(header.data)
end

function Base.done(header::SAMHeader, s)
    return done(header.data, s)
end

function Base.next(header::SAMHeader, s)
    return next(header.data, s)
end
