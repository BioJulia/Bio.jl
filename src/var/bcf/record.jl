# BCF Record
# ==========

type BCFRecord
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # fields
    sharedlen::UInt32
    indivlen::UInt32
end

function BCFRecord()
    return BCFRecord(UInt8[], 1:0, 0, 0)
end

function isfilled(record::BCFRecord)
    return !isempty(record.filled)
end

function datarange(record::BCFRecord)
    return record.filled
end

function checkfilled(record::BCFRecord)
    if !isfilled(record)
        throw(ArgumentError("unfilled BCF record"))
    end
end

function Base.:(==)(record1::BCFRecord, record2::BCFRecord)
    if isfilled(record1) == isfilled(record2) == true
        r1 = datarange(record1)
        r2 = datarange(record2)
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function Base.copy(record::BCFRecord)
    return BCFRecord(record.data[record.filled], record.filled, record.sharedlen, record.indivlen)
end

function BCFRecord(base::BCFRecord;
                   chromosome=nothing, position=nothing, quality=nothing,
                   identifier=nothing, reference=nothing, alternate=nothing,
                   filter=nothing, information=nothing, genotype=nothing)
    checkfilled(base)
    data = base.data[1:24]

    if chromosome != nothing
        store!(data, 0, convert(Int32, chromosome - 1))
    end
    if position != nothing
        store!(data, 4, convert(Int32, position - 1))
    end
    if quality != nothing
        store!(data, 12, convert(Float32, quality))
    end

    offset = boffset = 24
    if identifier == nothing
        offset, boffset = copyvec!(data, offset, base.data, boffset)
    else
        offset = storestr!(data, offset, string(identifier))
        boffset = skipvec(base.data, boffset)
    end

    if reference == nothing
        offset, boffset = copyvec!(data, offset, base.data, boffset)
    else
        reference = string(reference)
        offset = storestr!(data, offset, reference)
        if sizeof(reference) != rlen(base)
            store!(data, 8, convert(Int32, sizeof(reference)))
        end
        boffset = skipvec(base.data, boffset)
    end

    n = n_allele(base)
    if alternate == nothing
        for _ in 1:n-1
            offset, boffset = copyvec!(data, offset, base.data, boffset)
        end
    else
        if !isa(alternate, Vector)
            alternate = [alternate]
        end
        for x in alternate
            offset = storestr!(data, offset, string(x))
        end
        for _ in 1:n-1
            boffset = skipvec(base.data, boffset)
        end
        if length(alternate) != n
            n_allele_info = load(UInt32, data, 16)[1]
            n_allele_info = (n_allele_info & 0x0000ffff) | (UInt32(length(alternate) + 1) << 16)
            store!(data, 16, n_allele_info)
        end
    end

    if filter == nothing
        offset, boffset = copyvec!(data, offset, base.data, boffset)
    else
        if !isa(filter, Vector)
            filter = [filter]
        end
        offset = storevec!(data, offset, convert(Vector{Int8}, filter - 1))
        boffset = skipvec(base.data, boffset)
    end

    n = n_info(base)
    if information == nothing
        for _ in 1:n
            # copy key and value(s)
            offset, boffset = copyvec!(data, offset, base.data, boffset)
            offset, boffset = copyvec!(data, offset, base.data, boffset)
        end
    else
        if !isa(information, Associative)
            throw(ArgumentError("info must be an associative object"))
        end
        keyvec = Vector{Int8}(1)
        for (key, val) in information
            if !isa(key, Integer)
                throw(ArgumentError("info key must be an integer"))
            end
            keyvec[1] = key
            offset = storevec!(data, offset, keyvec)
            if isa(val, String)
                offset = storestr!(data, offset, val)
            elseif !isa(val, Vector)
                offset = storevec!(data, offset, [val])
            else
                offset = storevec!(data, offset, val)
            end
        end
        if length(information) != n
            if length(information) > typemax(UInt16)
                throw(ArgumentError("too many info fields"))
            end
            n_allele_info = load(UInt32, data, 16)[1]
            n_allele_info = (n_allele_info & 0xffff0000) | UInt32(length(information))
            store!(data, 16, n_allele_info)
        end
        for _ in 1:n
            # skip key and value(s)
            boffset = skipvec(base.data, boffset)
            boffset = skipvec(base.data, boffset)
        end
    end
    sharedlen = offset

    n = n_format(base)
    if genotype == nothing
        N = n_sample(base)
        for _ in 1:n
            # copy key
            offset, boffset = copyvec!(data, offset, base.data, boffset)
            # copy value(s)
            head, boffset = loadvechead(base.data, boffset)
            offset = storevechead!(data, offset, head)
            for _ in 1:N
                offset, boffset = copyvecbody!(data, offset, base.data, boffset, head)
            end
        end
    else
        error("modifying genotype is yet supported")
    end

    return BCFRecord(data, 1:endof(data), sharedlen, offset - sharedlen)
end

function Base.show(io::IO, record::BCFRecord)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   chromosome: ", chromosome(record))
        println(io, "     position: ", leftposition(record))
        println(io, "   identifier: ", identifier(record))
        println(io, "    reference: ", reference(record))
        println(io, "    alternate: ", join(alternate(record), ' '))
        println(io, "      quality: ", quality(record))
        println(io, "       filter: ", join(filter_(record), ' '))
          print(io, "  information: ", information(record))
    else
        print(io, " <not filled>")
    end
end

function chromosome(rec::BCFRecord)
    checkfilled(rec)
    return load(Int32, rec.data, 0)[1] % Int + 1
end

function leftposition(rec::BCFRecord)
    checkfilled(rec)
    return load(Int32, rec.data, 4)[1] % Int + 1
end

function rlen(rec::BCFRecord)
    checkfilled(rec)
    return load(Int32, rec.data, 8)[1] % Int
end

function quality(rec::BCFRecord)
    checkfilled(rec)
    # 0x7F800001 is a missing value.
    return load(Float32, rec.data, 12)[1]
end

function n_allele(rec::BCFRecord)
    checkfilled(rec)
    return (load(Int32, rec.data, 16)[1] >> 16) % Int
end

function n_info(rec::BCFRecord)
    checkfilled(rec)
    return (load(Int32, rec.data, 16)[1] & 0x0000ffff) % Int
end

function n_format(rec::BCFRecord)
    checkfilled(rec)
    return (load(UInt32, rec.data, 20)[1] >> 24) % Int
end

function n_sample(rec::BCFRecord)
    checkfilled(rec)
    return (load(UInt32, rec.data, 20)[1] & 0x000000ff) % Int
end

function identifier(rec::BCFRecord)
    checkfilled(rec)
    offset = 24
    return loadstr(rec.data, offset)[1]
end

function reference(rec::BCFRecord)
    checkfilled(rec)
    # skip ID
    offset = 24
    len, offset = loadveclen(rec.data, offset)
    # load REF
    return loadstr(rec.data, offset + len)[1]
end

function alternate(rec::BCFRecord)
    checkfilled(rec)
    # skip ID and REF
    offset = 24
    len, offset = loadveclen(rec.data, offset)
    len, offset = loadveclen(rec.data, offset + len)
    # load ALTs
    N = n_allele(rec) - 1
    alt = Vector{String}(N)
    for n in 1:N
        str, offset = loadstr(rec.data, offset + len)
        alt[n] = str
        len = sizeof(str)
    end
    return alt
end

function filter_(rec::BCFRecord)
    checkfilled(rec)
    # skip ID, REF and ALTs
    offset = 24
    len = 0
    for _ in 1:n_allele(rec)+1
        len, offset = loadveclen(rec.data, offset + len)
    end
    # load FILTER
    return loadvec(rec.data, offset + len)[1] .+ 1
end

function information(rec::BCFRecord; simplify::Bool=true)
    checkfilled(rec)
    # skip ID, REF, ALTs and FILTER
    offset::Int = 24
    len = 0
    for _ in 1:n_allele(rec)+2
        len, offset = loadveclen(rec.data, offset + len)
    end
    offset += len
    # load INFO
    ret = Vector{Tuple{Int,Any}}(n_info(rec))
    for i in 1:endof(ret)
        key, offset = loadvec(rec.data, offset)
        @assert length(key) == 1
        val, offset = loadvec(rec.data, offset)
        if simplify
            if isa(val, Vector) && length(val) == 1
                val = val[1]
            elseif isa(val, Vector{Void}) && isempty(val)
                val = nothing
            end
        end
        ret[i] = (key[1], val)
    end
    return ret
end

function information(rec::BCFRecord, key::Integer; simplify::Bool=true)
    checkfilled(rec)
    # skip ID, REF, ALTs and FILTER
    offset::Int = 24
    len = 0
    for _ in 1:n_allele(rec)+2
        len, offset = loadveclen(rec.data, offset + len)
    end
    offset += len
    # load INFO
    for _ in 1:n_info(rec)
        k, offset = loadvec(rec.data, offset)
        @assert length(k) == 1
        if k[1] == key
            val, offset = loadvec(rec.data, offset)
            if simplify
                if isa(val, Vector) && length(val) == 1
                    val = val[1]
                elseif isa(val, Vector{Void}) && isempty(val)
                    val = nothing
                end
            end
            return val
        else
            offset = skipvec(rec.data, offset)
        end
    end
    throw(KeyError(key))
end

function genotype(rec::BCFRecord)
    checkfilled(rec)
    offset::Int = rec.sharedlen
    N = n_sample(rec)
    ret = Tuple{Int,Any}[]
    for j in 1:n_format(rec)
        key, offset = loadvec(rec.data, offset)
        @assert length(key) == 1
        head, offset = loadvechead(rec.data, offset)
        vals = Vector{Any}[]
        for n in 1:N
            val, offset = loadvecbody(rec.data, offset, head)
            push!(vals, val)
        end
        push!(ret, (key[1], vals))
    end
    return ret
end

function genotype(rec::BCFRecord, index::Integer)
    return [(k, geno[index]) for (k, geno) in genotype(rec)]
end

function genotype(rec::BCFRecord, index::Integer, key::Integer)
    checkfilled(rec)
    N = n_sample(rec)
    offset::Int = rec.sharedlen
    for j in 1:n_format(rec)
        k, offset = loadvec(rec.data, offset)
        @assert length(k) == 1
        head, offset = loadvechead(rec.data, offset)
        for n in 1:N
            if k[1] == key && n == index
                return loadvecbody(rec.data, offset, head)[1]
            else
                offset = skipvecbody(rec.data, offset, head)
            end
        end
    end
    throw(KeyError(key))
end

function genotype{T<:Integer}(rec::BCFRecord, indexes::AbstractVector{T}, key::Integer)
    checkfilled(rec)
    N = n_sample(rec)
    offset::Int = rec.sharedlen
    for j in 1:n_format(rec)
        k, offset = loadvec(rec.data, offset)
        @assert length(k) == 1
        head, offset = loadvechead(rec.data, offset)
        if k[1] == key
            vals = Vector{Any}[]
            for n in 1:N
                val, offset = loadvecbody(rec.data, offset, head)
                push!(vals, val)
            end
            return vals[indexes]
        else
            for n in 1:N
                offset = skipvecbody(rec.data, offset, head)
            end
        end
    end
    throw(KeyError(key))
end

function genotype(rec::BCFRecord, ::Colon, key::Integer)
    return genotype(rec, 1:n_sample(rec), key)
end

function gt(x::Int8)
    allele = (x >> 1) - 1
    phased = x & 1 != 0
    return allele, phased
end

# Resolve type encoding.
function bcftype(b::UInt8)
    b &= 0x0f
    if b == 0x00
        # this is not defined in the BCF specs v2.2 but used in htslib
        return :void
    elseif b == 0x01
        return :int8
    elseif b == 0x02
        return :int16
    elseif b == 0x03
        return :int32
    elseif b == 0x05
        return :float32
    elseif b == 0x07
        return :character
    else
        return :unused
    end
end

function bcftypetag(typ::Symbol, len::Int)
    b = 0x00
    if typ == :void
        b |= 0x00
        @assert len == 0
    elseif typ == :int8
        b |= 0x01
    elseif typ == :int16
        b |= 0x02
    elseif typ == :int32
        b |= 0x03
    elseif typ == :float32
        b |= 0x05
    elseif typ == :character
        b |= 0x07
    else
        throw(ArgumentError("invalid type name: $(typ)"))
    end

    if len < 0
        throw(ArgumentError("length must be non-negative"))
    elseif len < 15
        b |= (len % UInt8) << 4
    else
        b |= (15 % UInt8)  << 4
    end
    return b
end

function bcftypesize(typ::Symbol)
    return typ == :int8      ? 1 :
           typ == :int16     ? 2 :
           typ == :int32     ? 4 :
           typ == :float32   ? 4 :
           typ == :character ? 1 :
           typ == :void      ? 0 :
           error("size unknown")
end

function load{T}(::Type{T}, data::Vector{UInt8}, offset::Int)
    if offset + sizeof(T) > length(data)
        throw(BoundsError())
    end
    return unsafe_load(convert(Ptr{T}, pointer(data, offset + 1))), offset + sizeof(T)
end

# Load a string; this is basically same as `loadvec` but type-stable.
function loadstr(data::Vector{UInt8}, offset::Int)
    (typ, len), offset = loadvechead(data, offset)
    @assert typ == :character
    return String(data[offset+1:offset+len]), offset
end

function loadvec(data::Vector{UInt8}, offset::Int)
    head, offset = loadvechead(data, offset)
    return loadvecbody(data, offset, head)
end

function loadveclen(data::Vector{UInt8}, offset::Int)
    (_, len), offset = loadvechead(data, offset)
    return len, offset
end

function loadvechead(data::Vector{UInt8}, offset::Int)
    b = data[offset+1]
    typ = bcftype(b)
    len::Int = b >> 4
    if len < 15
        return (typ, len), offset + 1
    end
    lenvec, offset = loadvec(data, offset + 1)
    @assert length(lenvec) == 1
    return (typ, Int(lenvec[1])), offset
end

function loadvecbody(data::Vector{UInt8}, offset::Int, head::Tuple{Symbol,Int})
    t, len = head
    T = t == :int8 ? Int8 :
        t == :int16 ? Int16 :
        t == :int32 ? Int32 :
        t == :float32 ? Float32 :
        t == :character ? UInt8 : Void
    ret = Vector{T}(len)
    for i in 1:len
        ret[i], offset = load(T, data, offset)
    end
    if t == :character
        return String(ret), offset
    else
        return ret, offset
    end
end

function skipvec(data::Vector{UInt8}, offset::Int)
    (t, len), offset = loadvechead(data, offset)
    return offset + bcftypesize(t) * len
end

function skipvecbody(data::Vector{UInt8}, offset::Int, head::Tuple{Symbol,Int})
    t, len = head
    return offset + bcftypesize(t) * len
end

function store!{T}(data::Vector{UInt8}, offset::Int, val::T)
    if offset + sizeof(T) > length(data)
        throw(BoundsError())
    end
    unsafe_store!(convert(Ptr{T}, pointer(data, offset + 1)), val)
    return offset + sizeof(T)
end

function storestr!(data::Vector{UInt8}, offset::Int, str::String)
    if !isascii(str)
        throw(ArgumentError("string must be ASCII"))
    end
    len = sizeof(str)
    offset = storevechead!(data, offset, (:character, len))
    if length(data) < offset + len
        resize!(data, offset + len)
    end
    memcpy(pointer(data, offset + 1), pointer(str), len)
    return offset + len
end

function storevec!{T}(data::Vector{UInt8}, offset::Int, vec::Vector{T})
    t = T == Int8 ? :int8 :
        T == Int16 ? :int16 :
        T == Int32 ? :int32 :
        T == Float32 ? :float32 :
        error("unsupported type: $(T)")
    len = length(vec)
    offset = storevechead!(data, offset, (t, len))
    return storevecbody!(data, offset, vec)
end

function storevechead!(data::Vector{UInt8}, offset::Int, head::Tuple{Symbol,Int})
    t, len = head
    if length(data) < offset + 1
        resize!(data, offset + 1)
    end
    data[offset+1] = bcftypetag(t, len)
    if len < 15
        return offset + 1
    elseif len < typemax(Int8)
        return storevec!(data, offset + 1, Int8[len])
    elseif len < typemax(Int16)
        return storevec!(data, offset + 1, Int16[len])
    elseif len < typemax(Int32)
        return storevec!(data, offset + 1, Int32[len])
    else
        error("too long vector")
    end
end

function storevecbody!{T}(data::Vector{UInt8}, offset::Int, vec::Vector{T})
    len = length(vec)
    n = len * sizeof(T)
    if length(data) < offset + n
        resize!(data, offset + n)
    end
    memcpy(pointer(data, offset + 1), pointer(vec), n)
    return offset + n
end

function copyvec!(dst::Vector{UInt8}, doffset::Int, src::Vector{UInt8}, soffset::Int)
    head, soffset = loadvechead(src, soffset)
    doffset = storevechead!(dst, doffset, head)
    return copyvecbody!(dst, doffset, src, soffset, head)
end

function copyvecbody!(dst::Vector{UInt8}, doffset::Int, src::Vector{UInt8}, soffset::Int, head::Tuple{Symbol,Int})
    t, len = head
    n = bcftypesize(t) * len
    if length(dst) < doffset + n
        resize!(dst, doffset + n)
    end
    copy!(dst, doffset + 1, src, soffset + 1, n)
    return doffset + n, soffset + n
end

function memcpy(dst::Ptr, src::Ptr, n::Int)
    ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Csize_t), dst, src, n)
    return nothing
end
