# Nucleotide Composition
# ======================

type NucleotideCounts{T <: Nucleotide}
    a::Uint
    c::Uint
    g::Uint
    t::Uint # also hold 'U' count when T == RNANucleotide
    n::Uint

    function NucleotideCounts()
        new(0, 0, 0, 0, 0)
    end
end


# Aliases
# -------

typealias DNANucleotideCounts NucleotideCounts{DNANucleotide}
typealias RNANucleotideCounts NucleotideCounts{RNANucleotide}


# Constructors
# ------------

# Count A, C, T/U, G respectively in a kmer stored in a Uint64
function count_a(x::Uint64)
    xinv = ~x
    return count_ones(((xinv >>> 1) & xinv) & 0x5555555555555555)
end
count_c(x::Uint64) = count_ones((((~x) >>> 1) & x) & 0x5555555555555555)
count_g(x::Uint64) = count_ones(((x >>> 1) & (~x)) & 0x5555555555555555)
count_t(x::Uint64) = count_ones((x    & (x >>> 1)) & 0x5555555555555555)

function NucleotideCounts{T}(seq::NucleotideSequence{T})
    dn, rn = divrem(seq.part.start - 1, 64)

    d = 2*dn
    r = 2*dn

    i = 1
    counts = NucleotideCounts{T}()

    # count leading unaligned bases
    for i in 1:r
        counts[seq[i]] += 1
        i += 1
    end
    if r > 0
        d += 1
    end

    # maybe just skip over blocks of Ns as I go?
    while i + 63 <= length(seq)
        # handle the all-N case
        if seq.ns.chunks[dn + 1] == 0xffffffffffffffff
            counts.n += 64
        else
            counts.a += count_a(seq.data[d + 1]) + count_a(seq.data[d + 2])
            counts.c += count_c(seq.data[d + 1]) + count_c(seq.data[d + 2])
            counts.g += count_g(seq.data[d + 1]) + count_g(seq.data[d + 2])
            counts.t += count_t(seq.data[d + 1]) + count_t(seq.data[d + 2])

            x = seq.ns.chunks[dn + 1]
            if x != 0
                for j in 1:64
                    if x & 0x01 != 0
                        counts.n += 1
                        counts[getnuc(T, seq.data, seq.part.start + i + j - 2)] -= 1
                    end

                    x >>= 1
                    if x == 0
                        break
                    end
                end
            end
        end

        dn += 1
        d += 2
        i += 64
    end

    # count trailing unaligned bases
    while i <= length(seq)
        counts[seq[i]] += 1
        i += 1
    end

    return counts
end

# Construct from K-mers
function NucleotideCounts{T,K}(seq::Kmer{T, K})
    x         = convert(Uint64, seq)
    counts    = NucleotideCounts{T}()
    counts.a += count_a(x) - 32 + K # Take leading zeros into account
    counts.c += count_c(x)
    counts.g += count_g(x)
    counts.t += count_t(x)
    return counts
end

# Basic Functions
# ---------------

getindex{T}(counts::NucleotideCounts{T}, nt::T) = getfield(counts, int(convert(Uint, nt) + 1))
setindex!{T}(counts::NucleotideCounts{T}, c::Integer, nt::T) = setfield!(counts, int(convert(Uint, nt) + 1), c)

# Pad strings so they are right justified when printed
function format_counts(xs)
    strings = String[string(x) for x in xs]
    len = maximum(map(length, strings))
    for i in 1:length(strings)
        strings[i] = string(repeat(" ", len - length(strings[i])), strings[i])
    end
    return strings
end


# Pretty printing of NucleotideCounts
function show(io::IO, counts::DNANucleotideCounts)
    count_strings = format_counts(
        [counts[DNA_A], counts[DNA_C], counts[DNA_G], counts[DNA_T], counts[DNA_N]])

    write(io,
        """
        DNANucleotideCounts:
          A => $(count_strings[1])
          C => $(count_strings[2])
          G => $(count_strings[3])
          T => $(count_strings[4])
          N => $(count_strings[5])
        """)
end


function show(io::IO, counts::RNANucleotideCounts)
    count_strings = format_counts(
        [counts[RNA_A], counts[RNA_C], counts[RNA_G], counts[RNA_U], counts[RNA_N]])

    write(io,
        """
        RNANucleotideCounts:
          A => $(count_strings[1])
          C => $(count_strings[2])
          G => $(count_strings[3])
          U => $(count_strings[4])
          N => $(count_strings[5])
        """)
end
