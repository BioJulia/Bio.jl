export readlsm

const _default_delims = [' ','\t']

"""
    readlsm([type,] io, delim=[' ', '\t'], comment='#')
    readlsm([type,] filepath, delim=[' ', '\t'], comment='#')

Reads a Labelled Square Matrix from an IO stream or file, returning a pair of
`(matrix, labels)`.

One may specify the element type of the resulting matrix; this defaults to
Float64.
"""
function readlsm{T}(::Type{T}, io::IO, delim=_default_delims, comment='#')
    matrix = Matrix{T}()
    labels = UTF8String[]
    nitems = 0
    row = 0
    for line in eachline(io)
        line = rstrip(line)
        if startswith(line, comment) || isempty(line)
            continue  # skip comments and empty lines
        end
        if row > nitems
            # Non-square error handled outside loop below
            break
        end
        cells = split(line, delim, keep=false)
        if isempty(labels)
            # header
            append!(labels, map(UTF8String, cells))
            nitems = length(labels)
            matrix = zeros(T, (nitems, nitems))
            continue
        end
        row += 1
        if cells[1] != labels[row]
            error("Cell label mismatch. (at row $row: $(cells[1]) != ",
                  "$(labels[row]))")
        end
        matrix[row, :] = map((s) -> parse(T, s), cells[2:end])
    end
    if row != nitems
        error("Non-square matrix: read $row rows, with $nitems columns.")
    end
    return (matrix, labels)
end

function readlsm(io::IO, delim=_default_delims)
    return readlsm(Float64, io, delim)
end

function readlsm{T}(::Type{T}, path::AbstractString, delim=_default_delims)
    open(path) do file
        return readlsm(T, file, delim)
    end
end

function readlsm(path::AbstractString, delim=_default_delims)
    return readlsm(Float64, path, delim)
end
