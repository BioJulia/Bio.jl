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

function Base.copy(rec::VCFRecord)
    return VCFRecord(
        copy(rec.data),
        rec.filled,
        rec.chrom,
        rec.pos,
        copy(rec.id),
        rec.ref,
        copy(rec.alt),
        rec.qual,
        copy(rec.filter),
        copy(rec.infokey),
        copy(rec.format),
        deepcopy(rec.genotype))
end

function Base.write(io::IO, rec::VCFRecord)
    checkfilled(rec)
    return write(io, rec.data)
end


# Accessor functions
# ------------------

"""
    chromosome(record::VCFRecord)::Nullable{String}

Get the chromosome name of `record`.
This function returns a null value if CHROM is missing.
"""
function chromosome(rec::VCFRecord)::Nullable{String}
    checkfilled(rec)
    if ismissing(rec.data, rec.chrom)
        return Nullable{String}()
    else
        return Nullable(String(rec.data[rec.chrom]))
    end
end

"""
    leftposition(record::VCFRecord)::Nullable{Int}

Get the reference position of `record`.
This function returns a null value if POS is missing.
"""
function Bio.leftposition(rec::VCFRecord)::Nullable{Int}
    checkfilled(rec)
    if ismissing(rec.data, rec.pos)
        return Nullable{Int}()
    else
        # TODO: no-copy accessor
        return Nullable(parse(Int, String(rec.data[rec.pos])))
    end
end

"""
    identifier(record::VCFRecord)::Vector{String}

Get the identifiers of `record`.
This function returns an empty vector if ID is missing.
"""
function identifier(rec::VCFRecord)::Vector{String}
    checkfilled(rec)
    if length(rec.id) == 1 && ismissing(rec.data, rec.id[1])
        return String[]
    else
        return [String(rec.data[r]) for r in rec.id]
    end
end

"""
    reference(record::VCFRecord)::Nullable{String}

Get the reference bases of `record`.
This function returns a null value if REF is missing.
"""
function reference(rec::VCFRecord)::Nullable{String}
    checkfilled(rec)
    if ismissing(rec.data, rec.ref)
        return Nullable{String}()
    else
        return Nullable(String(rec.data[rec.ref]))
    end
end

"""
    alternate(record::VCFRecord)::Vector{String}

Get the alternate bases of `record`.
This function returns an empty vector if ALT is missing.
"""
function alternate(rec::VCFRecord)::Vector{String}
    checkfilled(rec)
    if length(rec.alt) == 1 && ismissing(rec.data, rec.alt[1])
        return String[]
    else
        return [String(rec.data[r]) for r in rec.alt]
    end
end

"""
    quality(record::VCFRecord)::Nullable{Float64}

Get the quality score of `record`.
This function returns a null value if QUAL is missing.
"""
function quality(rec::VCFRecord)::Nullable{Float64}
    checkfilled(rec)
    if ismissing(rec.data, rec.qual)
        return Nullable{Float64}()
    else
        # TODO: no-copy parse
        return Nullable(parse(Float64, String(rec.data[rec.qual])))
    end
end

"""
    filter_(record::VCFRecord)::Vector{String}

Get the filter status of `record`.
This function returns an empty vector if FILTER is missing.
"""
function filter_(rec::VCFRecord)::Vector{String}
    checkfilled(rec)
    if length(rec.filter) == 1 && ismissing(rec.data, rec.filter[1])
        return String[]
    else
        return [String(rec.data[r]) for r in rec.filter]
    end
end

"""
    information(record::VCFRecord)::Vector{Pair{String,String}}

Get the additional information of `record`.
This function returns an empty vector if INFO is missing.
"""
function information(rec::VCFRecord)::Vector{Pair{String,String}}
    checkfilled(rec)
    ret = Pair{String,String}[]
    for (i, key) in enumerate(rec.infokey)
        val = infovalrange(rec, i)
        push!(ret, String(rec.data[key]) => String(rec.data[val]))
    end
    return ret
end

"""
    information(record::VCFRecord, key::String)::String

Get the additional information of `record` with `key`.
Keys without corresponding values return an empty string.
"""
function information(rec::VCFRecord, key::String)::String
    checkfilled(rec)
    # find key index
    i = 1
    while i ≤ endof(rec.infokey)
        if isequaldata(key, rec.data, rec.infokey[i])
            break
        end
        i += 1
    end
    if i > endof(rec.infokey)
        throw(KeyError(key))
    end
    val = infovalrange(rec, i)
    if isempty(val)
        return ""
    else
        return String(rec.data[val])
    end
end

"""
    infokeys(record::VCFRecord)::Vector{String}

Get the keys of the additional information of `record`.
This function returns an empty vector if INFO is missing.
"""
function infokeys(rec::VCFRecord)::Vector{String}
    checkfilled(rec)
    if length(rec.infokey) == 1 && ismissing(rec.data, rec.infokey[1])
        return String[]
    else
        return [String(rec.data[key]) for key in rec.infokey]
    end
end

# Returns the data range of the `i`-th value.
function infovalrange(rec::VCFRecord, i::Int)
    checkfilled(rec)
    data = rec.data
    key = rec.infokey[i]
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

Get the genotype format of `reocrd`.
This function returns an emtpy vector if FORMAT is missing.
"""
function format(rec::VCFRecord)::Vector{String}
    checkfilled(rec)
    if length(rec.format) == 1 && ismissing(rec.data, rec.format[1])
        return String[]
    else
        return [String(rec.data[r]) for r in rec.format]
    end
end

"""
    genotype(record::VCFRecord)::Vector{Vector{String}}

Get the genotypes of `record`.
"""
function genotype(rec::VCFRecord)
    checkfilled(rec)
    ret = Vector{String}[]
    for i in 1:endof(rec.genotype)
        push!(ret, genotype_impl(rec, i, 1:endof(rec.format)))
    end
    return ret
end

"""
    genotype(record::VCFRecord, index::Integer)::Vector{String}

Get the genotypes of the `index`-th individual in `record`.
This is effectively equivalent to `genotype(record)[index]` but more efficient.
"""
function genotype(rec::VCFRecord, index::Integer)
    checkfilled(rec)
    return genotype_impl(rec, index, 1:endof(rec.format))
end

"""
    genotype(record::VCFRecord, indexes, keys)

Get the genotypes in `record` that match `indexes` and `keys`.
`indexes` and `keys` can be either a scalar or a vector value.
"""
function genotype(rec::VCFRecord, index::Integer, key::String)::String
    checkfilled(rec)
    k = findgenokey(rec, key)
    if k == 0
        throw(KeyError(key))
    end
    return genotype_impl(rec, index, k)
end

function genotype(record::VCFRecord, index::Integer, keys::AbstractVector{String})::Vector{String}
    checkfilled(record)
    return [genotype(record, index, key) for key in keys]
end

function genotype{T<:Integer}(rec::VCFRecord, indexes::AbstractVector{T}, key::String)::Vector{String}
    checkfilled(rec)
    k = findgenokey(rec, key)
    if k == 0
        throw(KeyError(key))
    end
    return [genotype_impl(rec, i, k) for i in indexes]
end

function genotype{T<:Integer}(record::VCFRecord, indexes::AbstractVector{T}, keys::AbstractVector{String})::Vector{Vector{String}}
    checkfilled(rec)
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

function Base.show(io::IO, rec::VCFRecord)
    print(io, summary(rec), ':')
    if isfilled(rec)
        println(io)
        println(io, "   chromosome: ", get(chromosome(rec), "."))
        println(io, "     position: ", get(leftposition(rec), "."))
        println(io, "   identifier: ", let x = identifier(rec); isempty(x) ? "." : join(x, " "); end)
        println(io, "    reference: ", get(reference(rec), "."))
        println(io, "    alternate: ", let x = alternate(rec); isempty(x) ? "." : join(x, " "); end)
        println(io, "      quality: ", get(quality(rec), "."))
        println(io, "       filter: ", let x = filter_(rec); isempty(x) ? "." : join(x, " "); end)
        print(io, "  information: ")
        for (key, val) in information(rec)
            print(io, key)
            if !isempty(val)
                print(io, '=', val)
            end
            print(io, ' ')
        end
        println(io)
        println(io, "       format: ", let x = format(rec); isempty(x) ? "." : join(x, " "); end)
        print(io, "     genotype:")
        for i in 1:endof(rec.genotype)
            print(io, " [$(i)] ", let x = genotype(rec, i); isempty(x) ? "." : join(x, " "); end)
        end
    else
        print(io, " <not filled>")
    end
end

function ismissing(data::Vector{UInt8}, range::UnitRange{Int})
    return length(range) == 1 && data[first(range)] == UInt8('.')
end

# Check if `str == data[range]`
function isequaldata(str::String, data::Vector{UInt8}, range::UnitRange{Int})
    rlen = length(range)
    return rlen == sizeof(str) && memcmp(pointer(data, first(range)), pointer(str), rlen) == 0
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
end
