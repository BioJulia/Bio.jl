
# A Kmer is a sequence <= 32nt, without any 'N's, packed in a single 64 bit value.
#
# While NucleotideSequence is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.
#

# Note: while a Kmer isn't really an Unsigned, it's useful to declare it as such
# so Kmers can be used to intex into arrays.
bitstype 64 Kmer{T <: Nucleotide, K} <: Unsigned

typealias DNAKmer{K} Kmer{DNANucleotide, K}
typealias RNAKmer{K} Kmer{RNANucleotide, K}


# Conversion to/from Uint64
function convert(::Type{DNAKmer}, x::Uint64)
    return box(DNAKmer, unbox(Uint64, x))
end

function convert(::Type{RNAKmer}, x::Uint64)
    return box(RNAKmer, unbox(Uint64, x))
end

function convert(::Type{Uint64}, x::DNAKmer)
    return box(Uint64, unbox(DNAKmer, x))
end

function convert(::Type{Uint64}, x::RNAKmer)
    return box(Uint64, unbox(RNAKmer, x))
end


function show()

end


# TODO:
#   printing
#   conversion to/from string
#   conversion to/from NucleotideSequence
#   move eachkmer here
#


# TODO: Redo this. This needs to be as fast as possible. Benchmark, benchmark,
# benchmark.
# Call a function on every k-mer.
#function eachkmer{T}(f::Base.Callable, seq::NucleotideSequence{T}, k::Integer)
    #if k > 32
        #error("k must be ≤ 32 in eachkmer")
    #end
    #x = uint64(0)
    #mask = makemask(2 * k)
    #skip = k - 1
    #len = length(seq)
    #shift = 2 * (k - 1)
    #d, r = divrem(seq.part.start - 1, 64)

    ## manually iterate through N positions
    #ns_it = ns(seq)
    #ns_it_state = start(ns_it)
    #if done(ns_it, ns_it_state)
        #next_n_pos = length(seq) + 1
    #else
        #next_n_pos, ns_it_state = next(ns_it, ns_it_state)
    #end

    #i = 1
    #while i <= len
        ## skip over any kmer containing an N
        #if i == next_n_pos
            #if done(ns_it, ns_it_state)
                #next_n_pos = length(seq) + 1
            #else
                #next_n_pos, ns_it_state = next(ns_it, ns_it_state)
            #end
            #skip = k
        #end

        #x = (x >>> 2) | (((seq.data[d + 1] >>> r) & 0b11) << shift)

        #r += 2
        #if r >= 64
            #r = 0
            #d += 1
        #end

        #if skip == 0
            #f(i, x & mask)
            #skip = 1
        #else
            #skip -= 1
        #end
        #i += 1
    #end
#end



