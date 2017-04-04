
@generated function count_sites_bitpar{T<:Site,A}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})
    n = bitsof(A)
    if n == 2
        count_func = :count_bitpairs
        n_nucs = 32
    elseif n == 4
        count_func = :count_nibbles
        n_nucs = 16
    else
        error("n (= $n) ∉ (2, 4)")
    end

    quote
        if length(a) > length(b)
            return count_sites_bitpar(T, b, a)
        end
        @assert length(a) ≤ length(b)

        nexta = bitindex(a, 1)
        nextb = bitindex(b, 1)
        stopa = bitindex(a, endof(a) + 1)
        counts = start_count(T)
        spaces = 0

        # align reading position of `a.data` so that `offset(nexta) == 0`
        if nexta < stopa && offset(nexta) != 0
            #println("Aligning sequence A so it has no offset.")
            x = a.data[index(nexta)] >> offset(nexta)
            y = b.data[index(nextb)] >> offset(nextb)
            if offset(nextb) > offset(nexta)
                y |= b.data[index(nextb)+1] << (64 - offset(nextb))
            end
            k = 64 - offset(nexta)
            m = mask(k)
            #println("x masked: $(hex(x & m)), y masked: $(hex(y & m))")
            counts = update_count(counts, $count_func(T, x & m, y & m))
            #println("Count: ", counts)
            if correct_endspace(T)
                #println("N empty nibbles: ", nempty)
                nempty = $n_nucs - div(k, $n)
                counts = endspace_correction(nempty, counts)
                #println("Count, following empty correction: ", counts)
            end
            nexta += k
            nextb += k
            #println("Done initial alignment.")
        end
        @assert offset(nexta) == 0

        if offset(nextb) == 0  # data are aligned with each other
            #println("A and B are aligned with each other.")
            #println("Doing aligned loop.")
            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                #println("x: $(hex(x))")
                #println("y: $(hex(y))")
                counts = update_count(counts, $count_func(T, x, y))
                #println("Count: ", counts)
                nexta += 64
                nextb += 64
            end
            #println("Done aligned loop.")

            if nexta < stopa
                #println("Doing aligned tail.")
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                offs = stopa - nexta
                m = mask(offs)
                #println("x masked: $(hex(x & m))")
                #println("y masked: $(hex(y & m))")
                counts = update_count(counts, $count_func(T, x & m, y & m))
                #println("Count: ", counts)
                if correct_endspace(T)
                    nempty = $n_nucs - div(offs, $n)
                    #println("N empty nibbles: ", nempty)
                    counts = endspace_correction(nempty, counts)
                    #println("Count, following empty correction: ", counts)
                end
                #println("Done aligned tail.")
            end
        elseif nexta < stopa
            #println("A and B are not aligned.")
            #println("Doing unaligned loop.")

            y = b.data[index(nextb)]
            nextb += 64

            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                z = b.data[index(nextb)]
                y = y >> offset(nextb) | z << (64 - offset(nextb))

                #println("x: $(hex(x)), z: $(hex(z)), y: $(hex(y))")
                counts = update_count(counts, $count_func(T, x, y))
                #println("Count: ", counts)
                y = z
                nexta += 64
                nextb += 64
            end
            #println("Done unaligned loop.")

            if nexta < stopa
                #println("Doing unaligned tail.")
                x = a.data[index(nexta)]
                y = y >> offset(nextb)
                if 64 - offset(nextb) < stopa - nexta
                    y |= a.data[index(nextb)] << (64 - offset(nextb))
                end
                offs = stopa - nexta
                m = mask(offs)
                #println("x masked: $(hex(x & m)), y masked: $(hex(y & m))")
                counts = update_count(counts, $count_func(T, x & m, y & m))
                #println("Count: ", counts)
                if correct_endspace(T)
                    nempty = $n_nucs - div(offs, $n)
                    #println("N empty nibbles: ", nempty)
                    counts = endspace_correction(nempty, counts)
                    #println("Count, following empty correction: ", counts)
                end
                #println("Done unaligned tail.")
            end
        end
        #println("All done.")
        return counts
    end
end
