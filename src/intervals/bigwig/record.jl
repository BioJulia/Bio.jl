# BigBed Record
# =============

# type Record is defined in reader.jl.

"""
    BigWig.Record()

Create an unfilled bigWig record.
"""
function Record()
    return Record(0, 0, NaN32)
end

function Bio.isfilled(record::Record)
    return isdefined(record, :reader)
end

function Base.copy(record::Record)
    copy = Record(record.chromstart, record.chromend, record.value)
    copy.header = record.header
    if isdefined(record, :reader)
        copy.reader = record.reader
    end
    return copy
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "  chromosome: ", chrom(record))
        println(io, "       start: ", chromstart(record))
        println(io, "         end: ", chromend(record))
          print(io, "       value: ", value(record))
    else
        print(io, " <not filled>")
    end
end


# Accessor functions
# ------------------

"""
    chromid(record::Record)::UInt32

Get the chromosome ID of `record`.
"""
function chromid(record::Record)::UInt32
    checkfilled(record)
    return record.header.chrom_id
end

"""
    chrom(record::Record)::String

Get the chromosome name of `record`.
"""
function chrom(record::Record)::String
    checkfilled(record)
    return record.reader.chrom_names[chromid(record)]
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
