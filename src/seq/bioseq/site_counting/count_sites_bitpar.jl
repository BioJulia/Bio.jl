# count_sites_bitpar.jl
# =====================
#
# Counting sites in a bitparallel fashion.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

@inline bp_counter_type{S<:Site,A<:Alphabet}(::Type{S}, ::Type{A}) = Int
@inline bp_start_counter{S<:Site,A<:Alphabet}(::Type{S}, ::Type{A}) = zero(bp_counter_type(S, A))
@inline bp_update_counter(acc::Int, up::Int) = acc + up
@inline bp_correct_emptyspace{T<:Site,A<:Alphabet}(::Type{T}, ::Type{A}) = false
@inline bp_emptyspace_correction(nempty::Int, count::Int) = count - nempty

@generated function bitpar_counter{S<:Site,A<:NucleicAcidAlphabets}(::Type{S}, a::BioSequence{A}, b::BioSequence{A})
    n = bitsof(A)
    n_elems = div(64, n)

    quote
        # Ensure that `a` is always the shorter of the two sequences.
        if length(a) > length(b)
            return bitpar_counter(S, b, a)
        end
        @assert length(a) ≤ length(b)

        nexta = bitindex(a, 1)
        stopa = bitindex(a, endof(a) + 1)
        nextb = bitindex(b, 1)
        stopb = bitindex(b, endof(b) + 1)
        counts = bp_start_counter(S, A)

        println(A)
        println("Start indexes.")
        println("Index A: ", index(nexta), ", Offset A: ", offset(nexta))
        println("Index B: ", index(nextb), ", Offset B: ", offset(nextb))
        println("Stop Index A: ", index(stopa), ", Stop Offset A:", offset(stopa))
        println("Stop Index B: ", index(stopb), ", Stop Offset B:", offset(stopb))

        # The first thing we need to sort out is to correctly align the head of
        # sequence / subsequence `a`s data is aligned such that the offset of
        # `nexta` is essentially reduced to 0.
        # With sequence / subsequence `a` aligned, from there, we only need to
        # worry about the alignment of sequence / subsequence `b` with respect
        # to `a`.
        if nexta < stopa && offset(nexta) != 0

            println("Initial aligning, of a.data...")
            # Here we shift the first data chunks to the right so as the first
            # nucleotide of the seq/subseq is the first nibble / pair of bits.
            x = a.data[index(nexta)] >> offset(nexta)
            y = b.data[index(nextb)] >> offset(nextb)

            println("x: ", hex(x), "\ny: ", hex(y))

            # Here it is assumed that there is something to go and get from
            # the next chunk of `b`, yet that may not be true.
            # We know that if this is not true of `b`, then it is certainly not
            # true of `a`.
            println("stopb - nextb >= 64: ", stopb - nextb, " ", stopb - nextb > 64)
            if offset(nextb) > offset(nexta) && (stopb - nextb) > 64
                y |= b.data[index(nextb) + 1] << (64 - offset(nextb))
                println("Modified y: ", hex(y))
            end

            # Here we need to check something, we need to check if the
            # head of `a` we are currently aligning contains the entirity of
            # the sequence or subsequence.
            # Because if the end of seq/subseq `a` comes before the end of this
            # head integer, it's something we need to take into account of when
            # we mask x and y.
            #
            # In other words if stopa - nexta < 64, we know seq or subseq a's
            # data ends before the end of this data chunk, and so the mask
            # used needs to be defined to account for this: mask(stopa - nexta),
            # otherwise the mask simply needs to be mask(64 - offset(nexta)).
            #
            # This edge case was found and accounted for by Ben Ward @Ward9250.
            # As this maintainer for more information.
            println("stopa and nexta stuff:")
            println("stopa - nexta < 64: ", stopa - nexta < 64)
            println("stopa - nexta: ", stopa - nexta)
            println("stopa based mask: ", hex(mask(stopa - nexta)))
            println("64 - offset based mask: ", hex(mask(64 - offset(nexta))))

            k = ifelse(stopa - nexta < 64, stopa - nexta, 64 - offset(nexta))
            println("k: ", k)
            m = mask(k)
            println("mask used: ", hex(m))
            println("masked x: ", hex(x & m))
            println("masked y: ", hex(y & m))
            counts = bp_update_counter(counts, bp_chunk_count(S, A, x & m, y & m))
            if bp_correct_emptyspace(S, A)
                println("Correcting for emptyspace...")
                nempty = $n_elems - div(k, $n)
                println("nempty: ", nempty)
                counts = bp_emptyspace_correction(nempty, counts)
                println("counts: ", counts)
            end
            nexta += k
            nextb += k
            println("Index A: ", index(nexta), ", Offset A: ", offset(nexta))
            println("Index B: ", index(nextb), ", Offset B: ", offset(nextb))
            println("Finished Aligning of a.data...")
        end
        #@assert offset(nexta) == 0

        if offset(nextb) == 0  # data are aligned with each other
            println("Data are aligned...")
            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                println("x: ", hex(x))
                println("y: ", hex(y))
                counts = bp_update_counter(counts, bp_chunk_count(S, A, x, y))
                println("counts: ", counts)
                nexta += 64
                nextb += 64
            end

            if nexta < stopa
                println("Processing tail...")

                x = a.data[index(nexta)]
                println("x: ", hex(x))
                y = b.data[index(nextb)]
                println("y: ", hex(y))

                offs = stopa - nexta
                println("offs: ", offs)
                m = mask(offs)
                println("mask: ", hex(m))
                println("masked x: ", hex(x & m))
                println("masked y: ", hex(y & m))
                counts = bp_update_counter(counts, bp_chunk_count(S, A, x & m, y & m))
                println("counts: ", counts)
                if bp_correct_emptyspace(S, A)
                    println("Correcting for emptyspace...")
                    nempty = $n_elems - div(offs, $n)
                    println("nempty: ", nempty)
                    counts = bp_emptyspace_correction(nempty, counts)
                    println("counts: ", counts)
                end
            end
        elseif nexta < stopa
            println("Data are unaligned...")

            println("Index B: ", index(nextb), " Offset B: ", offset(nextb))
            y = b.data[index(nextb)]
            println("y: ", hex(y))
            nextb += 64
            # Note that here, updating `nextb` by 64, increases the chunk index,
            # but the `offset(nextb)` will remain the same.

            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                z = b.data[index(nextb)]
                y = y >> offset(nextb) | z << (64 - offset(nextb))
                println("x: ", hex(x))
                println("z: ", hex(z))
                println("y: ", hex(y))

                counts = bp_update_counter(counts, bp_chunk_count(S, A, x, y))
                println("counts: ", counts)
                y = z
                nexta += 64
                nextb += 64
            end

            if nexta < stopa
                println("Processing tail...")

                println("Index A: ", index(nexta), " Offset A: ", offset(nexta))
                println("Index B: ", index(nextb), " Offset B: ", offset(nextb))

                # Get the data chunk of `a`.
                x = a.data[index(nexta)]
                println("x: ", hex(x))
                y = y >> offset(nextb)
                println("y: ", hex(y))

                if 64 - offset(nextb) < stopa - nexta
                    #y |= a.data[index(nextb)] << (64 - offset(nextb))
                    y |= b.data[index(nextb)] << (64 - offset(nextb))
                end
                println("modified y: ", hex(y))

                offs = stopa - nexta
                println("offs: ", offs)
                m = mask(offs)
                println("mask: ", hex(m))
                println("masked x: ", hex(x & m))
                println("masked y: ", hex(y & m))
                counts = bp_update_counter(counts, bp_chunk_count(S, A, x & m, y & m))
                println("counts: ", counts)
                if bp_correct_emptyspace(S, A)
                    println("Correcting for emptyspace...")
                    nempty = $n_elems - div(offs, $n)
                    println("nempty: ", nempty)
                    counts = bp_emptyspace_correction(nempty, counts)
                    println("counts: ", counts)
                end
            end
        end
        return counts
    end
end
