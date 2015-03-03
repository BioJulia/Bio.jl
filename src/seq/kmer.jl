# K-mer
# =====


# A Kmer is a sequence <= 32nt, without any 'N's, packed in a single 64 bit value.
#
# While NucleotideSequence is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.

bitstype 64 Kmer{T<:Nucleotide, K}

typealias DNAKmer{K} Kmer{DNANucleotide, K}
typealias RNAKmer{K} Kmer{RNANucleotide, K}
typealias Codon RNAKmer{3}


# Conversion
# ----------

# Conversion to/from Uint64

convert{K}(::Type{DNAKmer{K}}, x::Uint64) = box(DNAKmer{K}, unbox(Uint64, x))
convert{K}(::Type{RNAKmer{K}}, x::Uint64) = box(RNAKmer{K}, unbox(Uint64, x))
convert(::Type{Uint64}, x::DNAKmer)       = box(Uint64, unbox(DNAKmer, x))
convert(::Type{Uint64}, x::RNAKmer)       = box(Uint64, unbox(RNAKmer, x))


# Convert to/from String

function convert{T}(::Type{Kmer{T}}, seq::String)
    k = length(seq)
    @assert k <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(Kmer{T, k}, seq)
end

function convert{T, K}(::Type{Kmer{T, K}}, seq::String)
    @assert length(seq) == K error("Cannot construct a $(K)-mer from a string of length $(length(seq))")

    x     = uint64(0)
    shift = 0
    for (i, c) in enumerate(seq)
        nt = convert(T, c)
        @assert nt != nnucleotide(T) error("A K-mer may not contain an N in its sequence")

        x |= convert(Uint64, nt) << shift
        shift += 2
    end

    return convert(Kmer{T, K}, x)
end

convert{T, K}(::Type{String}, seq::Kmer{T, K}) = convert(String, [convert(Char, x) for x in seq])


# Convert to/from NucleotideSequence
convert{T}(::Type{Kmer}, seq::NucleotideSequence{T}) = convert(Kmer{T}, seq)

function convert{T}(::Type{Kmer{T}}, seq::NucleotideSequence{T})
    k = length(seq)
    @assert k <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(Kmer{T, k}, seq)
end

function convert{T, K}(::Type{Kmer{T, K}}, seq::NucleotideSequence{T})
    @assert length(seq) == K error("Cannot construct a $(K)-mer from a NucleotideSequence of length $(length(seq))")

    x     = uint64(0)
    shift = 0
    for (i, nt) in enumerate(seq)
        if nt == nnucleotide(T)
            error("A K-mer may not contain an N in its sequence")
        end
        x |= convert(Uint64, nt) << shift
        shift += 2
    end
    return convert(Kmer{T, K}, x)
end

convert{T, K}(::Type{NucleotideSequence}, x::Kmer{T, K}) =  convert(NucleotideSequence{T}, x)

function convert{T, K}(::Type{NucleotideSequence{T}}, x::Kmer{T, K})
    ns = BitVector(K)
    fill!(ns, false)
    return NucleotideSequence{T}([convert(Uint64, x)], ns, 1:K)
end


# Constructors
# ------------

# From strings
dnakmer(seq::String) = convert(DNAKmer, seq)
rnakmer(seq::String) = convert(RNAKmer, seq)

function rnakmer(seq::String)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(RNAKmer{length(seq)}, seq)
end

# Constructors taking a sequence of nucleotides
function kmer{T <: Nucleotide}(nts::T...)
    K = length(nts)
    if K > 32
        error(string("Cannot construct a K-mer longer than 32nt."))
    end

    x = uint64(0)
    shift = 0
    for (i, nt) in enumerate(nts)
        if nt == nnucleotide(T)
            error("A Kmer may not contain an N in its sequence")
        end
        x |= convert(Uint64, nt) << shift
        shift += 2
    end
    return convert(Kmer{T, K}, x)
end

function dnakmer(seq::DNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(DNAKmer{length(seq)}, seq)
end


function dnakmer(seq::RNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(RNAKmer{length(seq)}, seq)
end


# Basic Functions
# ---------------

function =={T, K}(a::NucleotideSequence{T}, b::Kmer{T, K})
    if length(a) != K
        return false
    end

    for (u, v) in zip(a, b)
        if u != v
            return false
        end
    end

    return true
end


function =={T, K}(a::Kmer{T, K}, b::NucleotideSequence{T})
    return b == a
end

function getindex{T, K}(x::Kmer{T, K}, i::Integer)
    if i < 1 || i > K
        error(BoundsError())
    end
    return convert(T, (convert(Uint64, x) >>> (2*(i-1))) & 0b11)
end


function show{K}(io::IO, x::DNAKmer{K})
    write(io, "DNA $(K)-mer:\n ")
    for i in 1:K
        write(io, convert(Char, x[i]))
    end
end


function show{K}(io::IO, x::RNAKmer{K})
    write(io, "RNA $(K)-mer:\n ")
    for i in 1:K
        write(io, convert(Char, x[i]))
    end
end


isless{T, K}(x::Kmer{T, K}, y::Kmer{T, K}) = convert(Uint64, x) < convert(Uint64, y)

length{T, K}(x::Kmer{T, K}) = K

# Iterating over nucleotides
start(x::Kmer) = 1

function next{T, K}(x::Kmer{T, K}, i::Int)
    nt = convert(T, (convert(Uint64, x) >>> (2*(i-1))) & 0b11)
    return (nt, i + 1)
end

done{T, K}(x::Kmer{T, K}, i::Int) = i > K


# Other functions
# ---------------

reverse{T, K}(x::Kmer{T, K}) = convert(Kmer{T, K}, nucrev(convert(Uint64, x)) >>> (2 * (32 - K)))

function complement{T, K}(x::Kmer{T, K})
    return convert(Kmer{T, K},
        (~convert(Uint64, x)) & (0xffffffffffffffff >>> (2 * (32 - K))))
end

reverse_complement{T, K}(x::Kmer{T, K}) = complement(reverse(x))

mismatches{T, K}(x::Kmer{T, K}, y::Kmer{T, K}) = nucmismatches(convert(Uint64, x), convert(Uint64, y))


# A canonical k-mer is the numerical lesser of a k-mer and its reverse complement.
# This is useful in hashing/counting k-mers in data that is not strand specific,
# and thus observing k-mer is equivalent to observing its reverse complement.
function canonical{T, K}(x::Kmer{T, K})
    y = reverse_complement(x)
    return x < y ? x : y
end




# EachKmerIterator and EachKmerIteratorState
# ==========================================

# Iterate through every k-mer in a nucleotide sequence
immutable EachKmerIterator{T, K}
    seq::NucleotideSequence{T}
    nit::SequenceNIterator
    step::Int
end


immutable EachKmerIteratorState{T, K}
    i::Int
    x::Uint64
    next_n_pos::Int
    nit_state::Int
end


function eachkmer{T}(seq::NucleotideSequence{T}, k::Integer, step::Integer=1)
    if k < 0
        error("K must be ≥ 0 in EachKmer")
    elseif k > 32
        error("K must be ≤ 32 in EachKmer")
    end

    if step < 1
        error("step must be ≥ 1")
    end

    return EachKmerIterator{T, k}(seq, npositions(seq), step)
end


function nextkmer{T, K}(it::EachKmerIterator{T, K},
                        state::EachKmerIteratorState{T, K}, skip::Int)
    i = state.i + 1
    x = state.x
    next_n_pos = state.next_n_pos
    nit_state = state.nit_state

    shift = 2 * (K - 1)
    d, r = divrem(2 * (it.seq.part.start + i - 2), 64)
    while i <= length(it.seq)
        while next_n_pos < i
            if done(it.nit, nit_state)
                next_n_pos = length(it.seq) + 1
                break
            else
                next_n_pos, nit_state = next(it.nit, nit_state)
            end
        end

        if i - K + 1 <= next_n_pos <= i
            off = it.step * @compat ceil(Int, (K - skip) / it.step)
            if skip < K
                skip += off
            end
        end

        x = (x >>> 2) | (((it.seq.data[d + 1] >>> r) & 0b11) << shift)

        if skip == 0
            break
        end
        skip -= 1

        r += 2
        if r == 64
            r = 0
            d += 1
        end
        i += 1
    end

    return EachKmerIteratorState{T, K}(i, x, next_n_pos, nit_state)
end


function start{T, K}(it::EachKmerIterator{T, K})
    nit_state = start(it.nit)
    if done(it.nit, nit_state)
        next_n_pos = length(it.seq) + 1
    else
        next_n_pos, nit_state = next(it.nit, nit_state)
    end

    state = EachKmerIteratorState{T, K}(0, uint64(0), next_n_pos, nit_state)
    return nextkmer(it, state, K - 1)
end


function next{T, K}(it::EachKmerIterator{T, K},
                    state::EachKmerIteratorState{T, K})
    value = convert(Kmer{T, K}, state.x)
    next_state = nextkmer(it, state, it.step - 1)
    return (state.i - K + 1, value), next_state
end


function done{T, K}(it::EachKmerIterator{T, K},
                    state::EachKmerIteratorState{T, K})
    return state.i > length(it.seq)
end

# TODO: count_nucleotides
