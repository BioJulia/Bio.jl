
# A Kmer is a sequence <= 32nt, without any 'N's, packed in a single 64 bit value.
#
# While NucleotideSequence is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.
#

bitstype 64 Kmer{T <: Nucleotide, K}

typealias DNAKmer{K} Kmer{DNANucleotide, K}
typealias RNAKmer{K} Kmer{RNANucleotide, K}


# Conversion to/from Uint64
function convert{K}(::Type{DNAKmer{K}}, x::Uint64)
    return box(DNAKmer{K}, unbox(Uint64, x))
end

function convert{K}(::Type{RNAKmer{K}}, x::Uint64)
    return box(RNAKmer{K}, unbox(Uint64, x))
end

function convert(::Type{Uint64}, x::DNAKmer)
    return box(Uint64, unbox(DNAKmer, x))
end

function convert(::Type{Uint64}, x::RNAKmer)
    return box(Uint64, unbox(RNAKmer, x))
end


function convert{T}(::Type{Kmer{T}}, seq::String)
    k = length(seq)
    if k > 32
        error(string("Cannot construct a K-mer longer than 32nt."))
    end
    return convert(Kmer{T, k}, seq)
end


function convert{T, K}(::Type{Kmer{T, K}}, seq::String)
    if length(seq) != K
        error(string("Cannot construct a $(K)-mer from a string of length $(length(seq))"))
    end

    x = uint64(0)
    shift = 0
    for (i, c) in enumerate(seq)
        nt = convert(T, c)
        if nt == nnucleotide(T)
            error("A Kmer may not contain an N in its sequence")
        end
        x |= convert(Uint64, nt) << shift
        shift += 2
    end
    return convert(Kmer{T, K}, x)
end


function convert{T}(::Type{Kmer}, seq::NucleotideSequence{T})
    return convert(Kmer{T}, seq)
end


function convert{T}(::Type{Kmer{T}}, seq::NucleotideSequence{T})
    k = length(seq)
    if k > 32
        error(string("Cannot construct a K-mer longer than 32nt."))
    end
    return convert(Kmer{T, k}, seq)
end


function convert{T, K}(::Type{Kmer{T, K}}, seq::NucleotideSequence{T})
    if length(seq) != K
        error(string("Cannot construct a $(K)-mer from a string of length $(length(seq))"))
    end

    x = uint64(0)
    shift = 0
    for (i, nt) in enumerate(seq)
        if nt == nnucleotide(T)
            error("A Kmer may not contain an N in its sequence")
        end
        x |= convert(Uint64, nt) << shift
        shift += 2
    end
    return convert(Kmer{T, K}, x)
end


function convert{T, K}(::Type{NucleotideSequence}, x::Kmer{T, K})
    return convert(NucleotideSequence{T}, x)
end


function convert{T, K}(::Type{NucleotideSequence{T}}, x::Kmer{T, K})
    ns = BitVector(K)
    fill!(ns, false)
    return NucleotideSequence{T}([convert(Uint64, x)], ns, 1:K)
end


function convert{T, K}(::Type{String}, seq::Kmer{T, K})
    return convert(String, [convert(Char, x) for x in seq])
end


function dnakmer(seq::String)
    return convert(DNAKmer, seq)
end


function rnakmer(seq::String)
    return convert(RNAKmer, seq)
end


function dnakmer(seq::DNASequence)
    if length(seq) > 32
        error(string("Cannot construct a K-mer longer than 32nt."))
    end
    return convert(DNAKmer{length(seq)}, seq)
end


function dnakmer(seq::RNASequence)
    if length(seq) > 32
        error(string("Cannot construct a K-mer longer than 32nt."))
    end
    return convert(RNAKmer{length(seq)}, seq)
end


function rnakmer(seq::String)
    if length(seq) > 32
        error(string("Cannot construct a K-mer longer than 32nt."))
    end
    return convert(RNAKmer{length(seq)}, seq)
end


function getindex{T, K}(x::Kmer{T, K}, i::Integer)
    if i < 1 || i > K
        error(BoundsError())
    end
    convert(T, (convert(Uint64, x) >>> (2*(i-1))) & 0b11)
end


function show{T, K}(io::IO, x::Kmer{T, K})
    if T == DNANucleotide
        write(io, "dnakmer(\"")
    elseif T == RNANucleotide
        write(io, "rnakmer(\"")
    end

    for i in 1:K
        write(io, convert(Char, x[i]))
    end

    k = int(K)
    write(io, "\")  # ", T == DNANucleotide ? "DNA" : "RNA", " ", string(K), "-mer")
end


function isless{T, K}(x::Kmer{T, K}, y::Kmer{T, K})
    return convert(Uint64, x) < convert(Uint64, y)
end


function length{T, K}(x::Kmer{T, K})
    return K
end


# Iterate over nucleotides
function start(x::Kmer)
    return 1
end


function next{T, K}(x::Kmer{T, K}, i::Int)
    nt = convert(T, (convert(Uint64, x) >>> (2*(i-1))) & 0b11)
    return (nt, i + 1)
end


function done{T, K}(x::Kmer{T, K}, i::Int)
    return i > K
end


function reverse{T, K}(x::Kmer{T, K})
    return convert(Kmer{T, K}, nucrev(convert(Uint64, x)) >>> (2 * (32 - K)))
end


function complement{T, K}(x::Kmer{T, K})
    return convert(Kmer{T, K},
        (~convert(Uint64, x)) & (0xffffffffffffffff >>> (2 * (32 - K))))
end


function reverse_complement{T, K}(x::Kmer{T, K})
    return complement(reverse(x))
end


function mismatches{T, K}(x::Kmer{T, K}, y::Kmer{T, K})
    return nucmismatches(convert(Uint64, x), convert(Uint64, y))
end


# A canonical kmer is the numerical lesser of a k-mer and its reverse complement.
# This is useful in hashing/counting kmers in data that is not strand specific,
# and thus observing kmer is equivalent to observing its reverse complement.
function canonical{T, K}(x::Kmer{T, K})
    y = reverse_complement(x)
    return x < y ? x : y
end


# Iterate through every kmer in a nucleotide sequence
immutable EachKmerIterator{T, K}
    seq::NucleotideSequence{T}
    nit::SequenceNIterator
end


immutable EachKmerIteratorState{T, K}
    i::Int
    x::Uint64
    next_n_pos::Int
    nit_state::Int
end


function eachkmer{T}(seq::NucleotideSequence{T}, k::Integer)
    if k < 0
        error("K must be ≥ 0 in eachkmer")
    elseif k > 32
        error("K must be ≤ 32 in eachkmer")
    end

    return EachKmerIterator{T, k}(seq, npositions(seq))
end


function nextkmer{T, K}(it::EachKmerIterator{T, K},
                        state::EachKmerIteratorState{T, K}, skip)
    i = state.i + 1
    x = state.x
    next_n_pos = state.next_n_pos
    nit_state = state.nit_state

    shift = 2 * (K - 1)
    d, r = divrem(2 * (it.seq.part.start + i - 2), 64)
    while i <= length(it.seq)
        if i == next_n_pos
            if done(it.nit, nit_state)
                next_n_pos = length(it.seq) + 1
            else
                next_n_pos, nit_state = next(it.nit, nit_state)
            end
            skip = K
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
    #if state.i % 10000 == 0
        #println(STDERR, "next: ", state.i)
    #end
    value = convert(Kmer{T, K}, state.x)
    state = nextkmer(it, state, 0)
    return value, state
end


function done{T, K}(it::EachKmerIterator{T, K},
                    state::EachKmerIteratorState{T, K})
    return state.i > length(it.seq)
end

# TODO: count_nucleotides


