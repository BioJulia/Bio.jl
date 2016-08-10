module LabelledSquareMatrices

export readlsm

"""
Reads a Labelled Square Matrix

Returns a tuple of `(matrix, labels)`
"""
function readlsm{T}(::Type{T}, io::IO, delim='\t', comment='#')
    matrix = Matrix{T}()
    labels = UTF8String[]
    nitems = 0
    row = 0
    for line in eachline(io)
        line = chomp(line)
        if startswith(line, comment) || isempty(line)
            continue  # skip comments and empty lines
        end
        if row > nitems
            # Non-square error handled outside loop below
            break
        end
        cells = split(line, delim)
        if isempty(labels)
            # header
            append!(labels, map(UTF8String, cells[2:end]))
            nitems = length(labels)
            matrix = zeros(T, (nitems, nitems))
            continue
        end
        row += 1
        if cells[1] != labels[row]
            error("Cell label mismatch. (at row $row)")
        end
        matrix[row, :] = map((s) -> parse(T, s), cells[2:end])
    end
    if row != nitems
        error("Non-square matrix: read $row rows, with $nitems columns.")
    end
    return (matrix, labels)
end

function readlsm(io::IO, delim='\t')
    readlsm(Float64, io, delim)
end

function readlsm(path::AbstractString, delim='\t')
    open(path) do file
        lsm = readlsm(Float64, file, delim)
    end
end

end # module LabelledSquareMatrices
