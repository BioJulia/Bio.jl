
immutable BED <: Bio.IO.FileFormat end

export BED

function Base.open(filepath::AbstractString, mode::AbstractString, ::Type{BED};
                   n_fields::Integer=-1)
    io = open(filepath, mode)
    if mode[1] == 'r'
        return open(BufferedInputStream(io), BED)
    elseif mode[1] ∈ ('w', 'o')
        return BEDWriter(io, n_fields)
    end
    error("invalid open mode")
end

function Base.open(input::BufferedInputStream, ::Type{BED})
    return BEDReader(input)
end
