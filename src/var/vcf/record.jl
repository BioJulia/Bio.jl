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

function VCFRecord()
    return VCFRecord(
        # data and filled
        UInt8[], 1:0,
        # chrom-alt
        1:0, 1:0, UnitRange{Int}[], 1:0, UnitRange{Int}[],
        # qual-genotype
        1:0, UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[])
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

function Base.convert(::Type{VCFRecord}, str::String)
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

function chromosome(rec::VCFRecord)
    checkfilled(rec)
    if ismissing(rec.data, rec.chrom)
        return Nullable{String}()
    else
        return Nullable(String(rec.data[rec.chrom]))
    end
end

function leftposition(rec::VCFRecord)
    checkfilled(rec)
    if ismissing(rec.data, rec.pos)
        return Nullable{Int}()
    else
        return Nullable(parse(Int, String(rec.data[rec.pos])))
    end
end

function identifier(rec::VCFRecord)
    checkfilled(rec)
    if length(rec.id) == 1 && ismissing(rec.data, rec.id[1])
        return String[]
    else
        return [String(rec.data[r]) for r in rec.id]
    end
end

function reference(rec::VCFRecord)
    checkfilled(rec)
    if ismissing(rec.data, rec.ref)
        return Nullable{String}()
    else
        return Nullable(String(rec.data[rec.ref]))
    end
end

function alternate(rec::VCFRecord)
    checkfilled(rec)
    if length(rec.alt) == 1 && ismissing(rec.data, rec.alt[1])
        return String[]
    else
        return [String(rec.data[r]) for r in rec.alt]
    end
end

function quality(rec::VCFRecord)
    checkfilled(rec)
    if ismissing(rec.data, rec.qual)
        return Nullable{Float64}()
    else
        # TODO: no-copy parse
        return Nullable(parse(Float64, String(rec.data[rec.qual])))
    end
end

function filter_(rec::VCFRecord)
    checkfilled(rec)
    if length(rec.filter) == 1 && ismissing(rec.data, rec.filter[1])
        return String[]
    else
        return [String(rec.data[r]) for r in rec.filter]
    end
end

function information(rec::VCFRecord)
    checkfilled(rec)
    ret = Pair{String,String}[]
    for (i, key) in enumerate(rec.infokey)
        val = infovalrange(rec, i)
        push!(ret, String(rec.data[key]) => String(rec.data[val]))
    end
    return ret
end

function information(rec::VCFRecord, key::String)
    checkfilled(rec)
    i = findinfokey(rec, key)
    if i == 0
        throw(KeyError(key))
    end
    val = infovalrange(rec, i)
    if isempty(val)
        return ""
    else
        return String(rec.data[val])
    end
end

function findinfokey(rec::VCFRecord, key::String)
    for (i, infokey) in enumerate(rec.infokey)
        n = length(infokey)
        if n == sizeof(key) && memcmp(pointer(rec.data, first(infokey)), pointer(key), n) == 0
            return i
        end
    end
    return 0
end

function infokeys(rec::VCFRecord)
    checkfilled(rec)
    if length(rec.infokey) == 1 && ismissing(rec.data, rec.infokey[1])
        return String[]
    else
        return [String(rec.data[key]) for key in rec.infokey]
    end
end

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

function format(rec::VCFRecord)
    checkfilled(rec)
    if length(rec.format) == 1 && ismissing(rec.data, rec.format[1])
        return String[]
    else
        return [String(rec.data[r]) for r in rec.format]
    end
end

function genotype(rec::VCFRecord)
    checkfilled(rec)
    ret = Vector{String}[]
    for i in 1:endof(rec.genotype)
        push!(ret, genotype(rec, i))
    end
    return ret
end

function genotype(rec::VCFRecord, index::Integer)
    checkfilled(rec)
    ret = String[]
    indiv = rec.genotype[index]
    for j in 1:endof(rec.format)
        if j > endof(indiv)
            push!(ret, ".")  # missing
        else
            push!(ret, String(rec.data[indiv[j]]))
        end
    end
    return ret
end

function genotype(rec::VCFRecord, index::Integer, key::String)
    checkfilled(rec)
    j = findgeno(rec, key)
    if j == 0
        throw(KeyError(key))
    end
    indiv = rec.genotype[index]
    if endof(indiv) < j
        return "."
    else
        return String(rec.data[indiv[j]])
    end
end

function genotype{T<:Integer}(rec::VCFRecord, indexes::AbstractVector{T}, key::String)
    checkfilled(rec)
    j = findgeno(rec, key)
    if j == 0
        throw(KeyError(key))
    end
    ret = String[]
    for i in indexes
        indiv = rec.genotype[i]
        if endof(indiv) < j
            push!(ret, ".")  # missing
        else
            push!(ret, String(rec.data[indiv[j]]))
        end
    end
    return ret
end

function genotype(rec::VCFRecord, ::Colon, key::String)
    return genotype(rec, 1:endof(rec.genotype), key)
end

function findgeno(rec::VCFRecord, key::String)
    for (i, r) in enumerate(rec.format)
        n = length(r)
        if n == sizeof(key) && memcmp(pointer(rec.data, first(r)), pointer(key), n) == 0
            return i
        end
    end
    return 0
end

function ismissing(data::Vector{UInt8}, range::UnitRange{Int})
    return length(range) == 1 && data[first(range)] == UInt8('.')
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

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
end
