# BigBed Record
# =============

# type Record is defined in reader.jl.

function Record()
    return Record(Nullable{SectionHeader}(), 0, 0, NaN32)
end

function Bio.isfilled(record::Record)
    return !isnull(record.header)
end

function Base.copy(record::Record)
    return Record(record.header, record.chromstart, record.chromend, record.value)
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "  chromosome ID: ", chromid(record))
        println(io, "          start: ", chromstart(record))
        println(io, "            end: ", chromend(record))
          print(io, "          value: ", value(record))
    else
        print(io, " <not filled>")
    end
end


# Accessor functions
# ------------------

"""
    chromid(record::Record)::Int

Get the chromosome ID of `record`.
"""
function chromid(record::Record)::Int
    checkfilled(record)
    return get(record.header).chrom_id + 1
end

"""
    chromstart(record::Record)::Int

Get the start position of `record`.
"""
function chromstart(record::Record)::Int
    checkfilled(record)
    return record.chromstart + 1
end

"""
    chromend(record::Record)::Int

Get the end position of `record`.
"""
function chromend(record::Record)::Int
    checkfilled(record)
    return record.chromend % Int
end

"""
    value(record::Record)::Float32

Get the value of `record`.
"""
function value(record::Record)::Float32
    checkfilled(record)
    return record.value
end

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("bigWig record is not filled"))
    end
end
