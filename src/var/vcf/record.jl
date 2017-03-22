# VCF Record
# ==========

type Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # indexes
    chrom::UnitRange{Int}
    pos::UnitRange{Int}
    id::Vector{UnitRange{Int}}
    ref::UnitRange{Int}
    alt::Vector{UnitRange{Int}}
    qual::UnitRange{Int}
    filter::Vector{UnitRange{Int}}
    infokey::Vector{UnitRange{Int}}
    format::Vector{UnitRange{Int}}
    genotype::Vector{Vector{UnitRange{Int}}}
end

"""
    VCF.Record()

Create an unfilled VCF record.
"""
function Record()
    return Record(
        # data and filled
        UInt8[], 1:0,
        # chrom-alt
        1:0, 1:0, UnitRange{Int}[], 1:0, UnitRange{Int}[],
        # qual-genotype
        1:0, UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[])
end

"""
    VCF.Record(data::Vector{UInt8})

Create a VCF object from `data` containing a VCF record.
This function verifies the format and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(
        # data and filled
        data, 1:0,
        # chrom-alt
        1:0, 1:0, UnitRange{Int}[], 1:0, UnitRange{Int}[],
        # qual-genotype
        1:0, UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[])
    index!(record)
    return record
end

"""
    VCF.Record(str::AbstractString)

Create a VCF object from `str` containing a VCF record.
This function verifies the format and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return Record(convert(Vector{UInt8}, str))
end

function initialize!(record::Record)
    record.filled = 1:0
    record.chrom = 1:0
    record.pos = 1:0
    empty!(record.id)
    record.ref = 1:0
    empty!(record.alt)
    record.qual = 1:0
    empty!(record.filter)
    empty!(record.infokey)
    empty!(record.format)
    empty!(record.genotype)
    return record
end

function isfilled(record::Record)
    return !isempty(record.filled)
end

function datarange(record::Record)
    return record.filled
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = datarange(record1)
        r2 = datarange(record2)
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled VCF record"))
    end
end

function Record(base::Record;
                   chrom=nothing, pos=nothing, id=nothing,
                   ref=nothing, alt=nothing, qual=nothing,
                   filter=nothing, info=nothing, genotype=nothing)
    checkfilled(base)
    buf = IOBuffer()

    if chrom == nothing
        write(buf, base.data[base.chrom])
    else
        print(buf, string(chrom))
    end

    print(buf, '\t')
    if pos == nothing
        write(buf, base.data[base.pos])
    else
        print(buf, convert(Int, pos))
    end

    print(buf, '\t')
    if id == nothing
        if isempty(base.id)
            print(buf, '.')
        else
            for (i, r) in enumerate(base.id)
                if i != 1
                    print(buf, ';')
                end
                write(buf, base.data[r])
            end
        end
    else
        if !isa(id, Vector)
            id = [id]
        end
        if isempty(id)
            print(buf, '.')
        else
            for (i, x) in enumerate(id)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(x))
            end
        end
    end

    print(buf, '\t')
    if ref == nothing
        write(buf, base.data[base.ref])
    else
        print(buf, string(ref))
    end

    print(buf, '\t')
    if alt == nothing
        if isempty(base.alt)
            print(buf, '.')
        else
            for (i, r) in enumerate(base.alt)
                if i != 1
                    print(buf, ';')
                end
                write(buf, base.data[r])
            end
        end
    else
        if !isa(alt, Vector)
            alt = [alt]
        end
        if isempty(alt)
            print(buf, '.')
        else
            for (i, x) in enumerate(alt)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(x))
            end
        end
    end

    print(buf, '\t')
    if qual == nothing
        write(buf, base.data[base.qual])
    else
        print(buf, convert(Float64, qual))
    end

    print(buf, '\t')
    if filter == nothing
        if isempty(base.filter)
            print(buf, '.')
        else
            for (i, r) in enumerate(base.filter)
                if i != 1
                    print(buf, ';')
                end
                write(buf, base.data[r])
            end
        end
    else
        if !isa(filter, Vector)
            filter = [filter]
        end
        if isempty(filter)
            print(buf, '.')
        else
            for (i, x) in enumerate(filter)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(x))
            end
        end
    end

    print(buf, '\t')
    if info == nothing
        if isempty(base.infokey)
            print(buf, '.')
        else
            write(buf, base.data[first(base.infokey[1]):last(infovalrange(base, endof(base.infokey)))])
        end
    else
        if !isa(info, Associative)
            throw(ArgumentError("info must be an associative object"))
        elseif isempty(info)
            print(buf, '.')
        else
            for (i, (key, val)) in enumerate(info)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(key))
                if val != nothing
                    print(buf, '=', vcfformat(val))
                end
            end
        end
    end

    print(buf, '\t')
    if genotype == nothing
        if isempty(base.format)
            print(buf, '.')
        else
            write(buf, base.data[first(base.format[1]):last(base.format[end])])
        end
        if !isempty(base.genotype)
            for indiv in base.genotype
                print(buf, '\t')
                for (i, r) in enumerate(indiv)
                    if i != 1
                        print(buf, ':')
                    end
                    write(buf, base.data[r])
                end
            end
        end
    else
        if !isa(genotype, Vector)
            genotype = [genotype]
        end
        if isempty(genotype)
            print(buf, '.')
        else
            allkeys = String[]
            for indiv in genotype
                if !isa(indiv, Associative)
                    throw(ArgumentError("individual must be an associative"))
                end
                append!(allkeys, keys(indiv))
            end
            allkeys = sort!(unique(allkeys))
            if !isempty(allkeys)
                join(buf, allkeys, ':')
                for indiv in genotype
                    print(buf, '\t')
                    for (i, key) in enumerate(allkeys)
                        if i != 1
                            print(buf, ':')
                        end
                        print(buf, vcfformat(get(indiv, key, '.')))
                    end
                end
            end
        end
    end

    return Record(takebuf_array(buf))
end

function vcfformat(val)
    return string(val)
end

function vcfformat(val::Vector)
    return join(map(vcfformat, val), ',')
end

function Base.copy(record::Record)
    return Record(
        copy(record.data),
        record.filled,
        record.chrom,
        record.pos,
        copy(record.id),
        record.ref,
        copy(record.alt),
        record.qual,
        copy(record.filter),
        copy(record.infokey),
        copy(record.format),
        deepcopy(record.genotype))
end

function Base.write(io::IO, record::Record)
    checkfilled(record)
    return write(io, record.data)
end


# Accessor functions
# ------------------

"""
    chrom(record::Record)::String

Get the chromosome name of `record`.
"""
function chrom(record::Record)::String
    checkfilled(record)
    if ismissing(record, record.chrom)
        missingerror(:chrom)
    end
    return String(record.data[record.chrom])
end

function haschrom(record::Record)
    return isfilled(record) && !ismissing(record, record.chrom)
end

"""
    pos(record::Record)::Int

Get the reference position of `record`.
"""
function pos(record::Record)::Int
    checkfilled(record)
    if ismissing(record, record.pos)
        missingerror(:pos)
    end
    # TODO: no-copy accessor
    return parse(Int, String(record.data[record.pos]))
end

function haspos(record::Record)
    return isfilled(record) && !ismissing(record, record.pos)
end

"""
    id(record::Record)::Vector{String}

Get the identifiers of `record`.
"""
function id(record::Record)::Vector{String}
    checkfilled(record)
    if isempty(record.id)
        missingerror(:id)
    end
    return [String(record.data[r]) for r in record.id]
end

function hasid(record::Record)
    return isfilled(record) && !isempty(record.id)
end

"""
    ref(record::Record)::String

Get the reference bases of `record`.
"""
function ref(record::Record)::String
    checkfilled(record)
    if ismissing(record, record.ref)
        missingerror(:ref)
    end
    return String(record.data[record.ref])
end

function hasref(record::Record)
    return isfilled(record) && !ismissing(record, record.ref)
end

"""
    alt(record::Record)::Vector{String}

Get the alternate bases of `record`.
"""
function alt(record::Record)::Vector{String}
    checkfilled(record)
    if isempty(record.alt)
        missingerror(:alt)
    end
    return [String(record.data[r]) for r in record.alt]
end

function hasalt(record::Record)
    return isfilled(record) && !isempty(record.alt)
end

"""
    qual(record::Record)::Float64

Get the quality score of `record`.
"""
function qual(record::Record)::Float64
    checkfilled(record)
    if ismissing(record, record.qual)
        missingerror(:qual)
    end
    # TODO: no-copy parse
    return parse(Float64, String(record.data[record.qual]))
end

function hasqual(record::Record)
    return isfilled(record) && !ismissing(record, record.qual)
end

"""
    filter(record::Record)::Vector{String}

Get the filter status of `record`.
"""
function filter(record::Record)::Vector{String}
    checkfilled(record)
    if isempty(record.filter)
        missingerror(:filter)
    end
    return [String(record.data[r]) for r in record.filter]
end

function hasfilter(record::Record)
    return isfilled(record) && !isempty(record.filter)
end

"""
    info(record::Record)::Vector{Pair{String,String}}

Get the additional information of `record`.
"""
function info(record::Record)::Vector{Pair{String,String}}
    checkfilled(record)
    if isempty(record.infokey)
        missingerror(:info)
    end
    ret = Pair{String,String}[]
    for (i, key) in enumerate(record.infokey)
        val = infovalrange(record, i)
        push!(ret, String(record.data[key]) => String(record.data[val]))
    end
    return ret
end

function hasinfo(record::Record)
    return isfilled(record) && !isempty(record.infokey)
end

"""
    info(record::Record, key::String)::String

Get the additional information of `record` with `key`.
Keys without corresponding values return an empty string.
"""
function info(record::Record, key::String)::String
    checkfilled(record)
    i = findinfokey(record, key)
    if i == 0
        throw(KeyError(key))
    end
    val = infovalrange(record, i)
    if isempty(val)
        return ""
    else
        return String(record.data[val])
    end
end

function hasinfo(record::Record, key::String)
    return isfilled(record) && findinfokey(key) > 0
end

function findinfokey(record::Record, key::String)
    for i in 1:endof(record.infokey)
        if isequaldata(key, record.data, record.infokey[i])
            return i
        end
    end
    return 0
end

"""
    infokeys(record::Record)::Vector{String}

Get the keys of the additional information of `record`.
This function returns an empty vector when the INFO field is missing.
"""
function infokeys(record::Record)::Vector{String}
    checkfilled(record)
    return [String(record.data[key]) for key in record.infokey]
end

# Returns the data range of the `i`-th value.
function infovalrange(record::Record, i::Int)
    checkfilled(record)
    data = record.data
    key = record.infokey[i]
    if last(key) + 1 ≤ endof(data) && data[last(key)+1] == UInt8('=')
        endpos = search(data, ';', last(key)+1)
        if endpos == 0
            endpos = search(data, '\t', last(key)+1)
            @assert endpos != 0
        end
        return last(key)+2:endpos-1
    else
        return last(key)+1:last(key)
    end
end

"""
    format(record::Record)::Vector{String}

Get the genotype format of `record`.
"""
function format(record::Record)::Vector{String}
    checkfilled(record)
    if isempty(record.format)
        missingerror(:format)
    end
    return [String(record.data[r]) for r in record.format]
end

function hasformat(record::Record)
    return isfilled(record) && !isempty(record.format)
end

"""
    genotype(record::Record)::Vector{Vector{String}}

Get the genotypes of `record`.
"""
function genotype(record::Record)
    checkfilled(record)
    ret = Vector{String}[]
    for i in 1:endof(record.genotype)
        push!(ret, genotype_impl(record, i, 1:endof(record.format)))
    end
    return ret
end

"""
    genotype(record::Record, index::Integer)::Vector{String}

Get the genotypes of the `index`-th individual in `record`.
This is effectively equivalent to `genotype(record)[index]` but more efficient.
"""
function genotype(record::Record, index::Integer)
    checkfilled(record)
    return genotype_impl(record, index, 1:endof(record.format))
end

"""
    genotype(record::Record, indexes, keys)

Get the genotypes in `record` that match `indexes` and `keys`.
`indexes` and `keys` can be either a scalar or a vector value.
Trailing fields that are dropped are filled with `"."`.
"""
function genotype(record::Record, index::Integer, key::String)::String
    checkfilled(record)
    k = findgenokey(record, key)
    if k == 0
        throw(KeyError(key))
    end
    return genotype_impl(record, index, k)
end

function genotype(record::Record, index::Integer, keys::AbstractVector{String})::Vector{String}
    checkfilled(record)
    return [genotype(record, index, key) for key in keys]
end

function genotype{T<:Integer}(record::Record, indexes::AbstractVector{T}, key::String)::Vector{String}
    checkfilled(record)
    k = findgenokey(record, key)
    if k == 0
        throw(KeyError(key))
    end
    return [genotype_impl(record, i, k) for i in indexes]
end

function genotype{T<:Integer}(record::Record, indexes::AbstractVector{T}, keys::AbstractVector{String})::Vector{Vector{String}}
    checkfilled(record)
    ks = Vector{Int}(length(keys))
    for i in 1:endof(keys)
        key = keys[i]
        k = findgenokey(record, key)
        if k == 0
            throw(KeyError(key))
        end
        ks[i] = k
    end
    return [genotype_impl(record, i, ks) for i in indexes]
end

function genotype(record::Record, ::Colon, key::String)::Vector{String}
    return genotype(record, 1:endof(record.genotype), key)
end

function findgenokey(record::Record, key::String)
    return findfirst(r -> isequaldata(key, record.data, r), record.format)
end

function genotype_impl(record::Record, index::Int, keys::AbstractVector{Int})
    return [genotype_impl(record, index, k) for k in keys]
end

function genotype_impl(record::Record, index::Int, key::Int)
    geno = record.genotype[index]
    if key > endof(geno)  # dropped field
        return "."
    else
        return String(record.data[geno[key]])
    end
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   chromosome: ", haschrom(record) ? chrom(record) : "<missing>")
        println(io, "     position: ", haspos(record) ? pos(record) : "<missing>")
        println(io, "   identifier: ", hasid(record) ? join(id(record), " ") : "<missing>")
        println(io, "    reference: ", hasref(record) ? ref(record) : "<missing>")
        println(io, "    alternate: ", hasalt(record) ? join(alt(record), " ") : "<missing>")
        println(io, "      quality: ", hasqual(record) ? qual(record) : "<missing>")
        println(io, "       filter: ", hasfilter(record) ? join(filter(record), " ") : "<missing>")
          print(io, "  information: ")
        if hasinfo(record)
            for (key, val) in info(record)
                print(io, key)
                if !isempty(val)
                    print(io, '=', val)
                end
                print(io, ' ')
            end
        else
            print(io, "<missing>")
        end
        println(io)
          print(io, "       format: ", hasformat(record) ? join(format(record), " ") : "<missing>")
        if hasformat(record)
            println(io)
            print(io, "     genotype:")
            for i in 1:endof(record.genotype)
                print(io, " [$(i)] ", let x = genotype(record, i); isempty(x) ? "." : join(x, " "); end)
            end
        end
    else
        print(io, " <not filled>")
    end
end

function ismissing(record::Record, range::UnitRange{Int})
    return length(range) == 1 && record.data[first(range)] == UInt8('.')
end

# Check if `str == data[range]`
function isequaldata(str::String, data::Vector{UInt8}, range::UnitRange{Int})
    rlen = length(range)
    return rlen == sizeof(str) && memcmp(pointer(data, first(range)), pointer(str), rlen) == 0
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
end
