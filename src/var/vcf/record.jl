# VCF Record
# ==========

type VCFRecord
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
    VCFRecord()

Create an unfilled `VCFRecord` object.
"""
function VCFRecord()
    return VCFRecord(
        # data and filled
        UInt8[], 1:0,
        # chrom-alt
        1:0, 1:0, UnitRange{Int}[], 1:0, UnitRange{Int}[],
        # qual-genotype
        1:0, UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[])
end

"""
    VCFRecord(data::Vector{UInt8})

Create a `VCFRecord` object from `data` containing a VCF record.
This function verifies the format and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new `VCFRecord` object.
"""
function VCFRecord(data::Vector{UInt8})
    return convert(VCFRecord, data)
end

function Base.convert(::Type{VCFRecord}, data::Vector{UInt8})
    record = VCFRecord(
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
    VCFRecord(str::AbstractString)

Create a `VCFRecord` object from `str` containing a VCF record.
This function verifies the format and indexes fields for accessors.
"""
function VCFRecord(str::AbstractString)
    return convert(VCFRecord, str)
end

function Base.convert(::Type{VCFRecord}, str::AbstractString)
    return VCFRecord(convert(Vector{UInt8}, str))
end

function initialize!(record::VCFRecord)
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

function isfilled(record::VCFRecord)
    return !isempty(record.filled)
end

function datarange(record::VCFRecord)
    return record.filled
end

function Base.:(==)(record1::VCFRecord, record2::VCFRecord)
    if isfilled(record1) == isfilled(record2) == true
        r1 = datarange(record1)
        r2 = datarange(record2)
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function checkfilled(record::VCFRecord)
    if !isfilled(record)
        throw(ArgumentError("unfilled VCF record"))
    end
end

function VCFRecord(base::VCFRecord;
                   chromosome=nothing, position=nothing, identifier=nothing,
                   reference=nothing, alternate=nothing, quality=nothing,
                   filter=nothing, information=nothing, genotype=nothing)
    checkfilled(base)
    buf = IOBuffer()

    if chromosome == nothing
        write(buf, base.data[base.chrom])
    else
        print(buf, string(chromosome))
    end

    print(buf, '\t')
    if position == nothing
        write(buf, base.data[base.pos])
    else
        print(buf, convert(Int, position))
    end

    print(buf, '\t')
    if identifier == nothing
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
        if !isa(identifier, Vector)
            identifier = [identifier]
        end
        if isempty(identifier)
            print(buf, '.')
        else
            for (i, x) in enumerate(identifier)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(x))
            end
        end
    end

    print(buf, '\t')
    if reference == nothing
        write(buf, base.data[base.ref])
    else
        print(buf, string(reference))
    end

    print(buf, '\t')
    if alternate == nothing
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
        if !isa(alternate, Vector)
            alternate = [alternate]
        end
        if isempty(alternate)
            print(buf, '.')
        else
            for (i, x) in enumerate(alternate)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(x))
            end
        end
    end

    print(buf, '\t')
    if quality == nothing
        write(buf, base.data[base.qual])
    else
        print(buf, convert(Float64, quality))
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
    if information == nothing
        if isempty(base.infokey)
            print(buf, '.')
        else
            write(buf, base.data[first(base.infokey[1]):last(infovalrange(base, endof(base.infokey)))])
        end
    else
        if !isa(information, Associative)
            throw(ArgumentError("information must be an associative object"))
        elseif isempty(information)
            print(buf, '.')
        else
            for (i, (key, val)) in enumerate(information)
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

    return VCFRecord(takebuf_array(buf))
end

function vcfformat(val)
    return string(val)
end

function vcfformat(val::Vector)
    return join(map(vcfformat, val), ',')
end

function Base.copy(record::VCFRecord)
    return VCFRecord(
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

function Base.write(io::IO, record::VCFRecord)
    checkfilled(record)
    return write(io, record.data)
end


# Accessor functions
# ------------------

"""
    chromosome(record::VCFRecord)::String

Get the chromosome name of `record`.
"""
function chromosome(record::VCFRecord)::String
    checkfilled(record)
    if ismissing(record, record.chrom)
        missingerror(:chromosome)
    end
    return String(record.data[record.chrom])
end

function haschromosome(record::VCFRecord)
    return isfilled(record) && !ismissing(record, record.chrom)
end

"""
    leftposition(record::VCFRecord)::Int

Get the reference position of `record`.
"""
function Bio.leftposition(record::VCFRecord)::Int
    checkfilled(record)
    if ismissing(record, record.pos)
        missingerror(:leftposition)
    end
    # TODO: no-copy accessor
    return parse(Int, String(record.data[record.pos]))
end

function Bio.hasleftposition(record::VCFRecord)
    return isfilled(record) && !ismissing(record, record.pos)
end

"""
    identifier(record::VCFRecord)::Vector{String}

Get the identifiers of `record`.
"""
function identifier(record::VCFRecord)::Vector{String}
    checkfilled(record)
    if isempty(record.id)
        missingerror(:identifier)
    end
    return [String(record.data[r]) for r in record.id]
end

function hasidentifier(record::VCFRecord)
    return isfilled(record) && !isempty(record.id)
end

"""
    reference(record::VCFRecord)::String

Get the reference bases of `record`.
"""
function reference(record::VCFRecord)::String
    checkfilled(record)
    if ismissing(record, record.ref)
        missingerror(:reference)
    end
    return String(record.data[record.ref])
end

function hasreference(record::VCFRecord)
    return isfilled(record) && !ismissing(record, record.ref)
end

"""
    alternate(record::VCFRecord)::Vector{String}

Get the alternate bases of `record`.
"""
function alternate(record::VCFRecord)::Vector{String}
    checkfilled(record)
    if isempty(record.alt)
        missingerror(:alternate)
    end
    return [String(record.data[r]) for r in record.alt]
end

function hasalternate(record::VCFRecord)
    return isfilled(record) && !isempty(record.alt)
end

"""
    quality(record::VCFRecord)::Float64

Get the quality score of `record`.
"""
function quality(record::VCFRecord)::Float64
    checkfilled(record)
    if ismissing(record, record.qual)
        missingerror(:quality)
    end
    # TODO: no-copy parse
    return parse(Float64, String(record.data[record.qual]))
end

function hasquality(record::VCFRecord)
    return isfilled(record) && !ismissing(record, record.qual)
end

"""
    filter_(record::VCFRecord)::Vector{String}

Get the filter status of `record`.
"""
function filter_(record::VCFRecord)::Vector{String}
    checkfilled(record)
    if isempty(record.filter)
        missingerror(:filter_)
    end
    return [String(record.data[r]) for r in record.filter]
end

function hasfilter_(record::VCFRecord)
    return isfilled(record) && !isempty(record.filter)
end

"""
    information(record::VCFRecord)::Vector{Pair{String,String}}

Get the additional information of `record`.
"""
function information(record::VCFRecord)::Vector{Pair{String,String}}
    checkfilled(record)
    if isempty(record.infokey)
        missingerror(:information)
    end
    ret = Pair{String,String}[]
    for (i, key) in enumerate(record.infokey)
        val = infovalrange(record, i)
        push!(ret, String(record.data[key]) => String(record.data[val]))
    end
    return ret
end

function hasinformation(record::VCFRecord)
    return isfilled(record) && !isempty(record.infokey)
end

"""
    information(record::VCFRecord, key::String)::String

Get the additional information of `record` with `key`.
Keys without corresponding values return an empty string.
"""
function information(record::VCFRecord, key::String)::String
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

function hasinformation(record::VCFRecord, key::String)
    return isfilled(record) && findinfokey(key) > 0
end

function findinfokey(record::VCFRecord, key::String)
    for i in 1:endof(record.infokey)
        if isequaldata(key, record.data, record.infokey[i])
            return i
        end
    end
    return 0
end

"""
    infokeys(record::VCFRecord)::Vector{String}

Get the keys of the additional information of `record`.
This function returns an empty vector when the INFO field is missing.
"""
function infokeys(record::VCFRecord)::Vector{String}
    checkfilled(record)
    return [String(record.data[key]) for key in record.infokey]
end

# Returns the data range of the `i`-th value.
function infovalrange(record::VCFRecord, i::Int)
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
    format(record::VCFRecord)::Vector{String}

Get the genotype format of `record`.
"""
function format(record::VCFRecord)::Vector{String}
    checkfilled(record)
    if isempty(record.format)
        missingerror(:format)
    end
    return [String(record.data[r]) for r in record.format]
end

function hasformat(record::VCFRecord)
    return isfilled(record) && !isempty(record.format)
end

"""
    genotype(record::VCFRecord)::Vector{Vector{String}}

Get the genotypes of `record`.
"""
function genotype(record::VCFRecord)
    checkfilled(record)
    ret = Vector{String}[]
    for i in 1:endof(record.genotype)
        push!(ret, genotype_impl(record, i, 1:endof(record.format)))
    end
    return ret
end

"""
    genotype(record::VCFRecord, index::Integer)::Vector{String}

Get the genotypes of the `index`-th individual in `record`.
This is effectively equivalent to `genotype(record)[index]` but more efficient.
"""
function genotype(record::VCFRecord, index::Integer)
    checkfilled(record)
    return genotype_impl(record, index, 1:endof(record.format))
end

"""
    genotype(record::VCFRecord, indexes, keys)

Get the genotypes in `record` that match `indexes` and `keys`.
`indexes` and `keys` can be either a scalar or a vector value.
Trailing fields that are dropped are filled with `"."`.
"""
function genotype(record::VCFRecord, index::Integer, key::String)::String
    checkfilled(record)
    k = findgenokey(record, key)
    if k == 0
        throw(KeyError(key))
    end
    return genotype_impl(record, index, k)
end

function genotype(record::VCFRecord, index::Integer, keys::AbstractVector{String})::Vector{String}
    checkfilled(record)
    return [genotype(record, index, key) for key in keys]
end

function genotype{T<:Integer}(record::VCFRecord, indexes::AbstractVector{T}, key::String)::Vector{String}
    checkfilled(record)
    k = findgenokey(record, key)
    if k == 0
        throw(KeyError(key))
    end
    return [genotype_impl(record, i, k) for i in indexes]
end

function genotype{T<:Integer}(record::VCFRecord, indexes::AbstractVector{T}, keys::AbstractVector{String})::Vector{Vector{String}}
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

function genotype(record::VCFRecord, ::Colon, key::String)::Vector{String}
    return genotype(record, 1:endof(record.genotype), key)
end

function findgenokey(record::VCFRecord, key::String)
    return findfirst(r -> isequaldata(key, record.data, r), record.format)
end

function genotype_impl(record::VCFRecord, index::Int, keys::AbstractVector{Int})
    return [genotype_impl(record, index, k) for k in keys]
end

function genotype_impl(record::VCFRecord, index::Int, key::Int)
    geno = record.genotype[index]
    if key > endof(geno)  # dropped field
        return "."
    else
        return String(record.data[geno[key]])
    end
end

function Base.show(io::IO, record::VCFRecord)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   chromosome: ", haschromosome(record) ? chromosome(record) : "<missing>")
        println(io, "     position: ", hasleftposition(record) ? leftposition(record) : "<missing>")
        println(io, "   identifier: ", hasidentifier(record) ? join(identifier(record), " ") : "<missing>")
        println(io, "    reference: ", hasreference(record) ? reference(record) : "<missing>")
        println(io, "    alternate: ", hasalternate(record) ? join(alternate(record), " ") : "<missing>")
        println(io, "      quality: ", hasquality(record) ? quality(record) : "<missing>")
        println(io, "       filter: ", hasfilter_(record) ? join(filter_(record), " ") : "<missing>")
          print(io, "  information: ")
        if hasinformation(record)
            for (key, val) in information(record)
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

function ismissing(record::VCFRecord, range::UnitRange{Int})
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
