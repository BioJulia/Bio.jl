BC = BitparCount

@inline correct_endspace{T<:Site,A<:Alphabet}(::Type{T}, ::Type{A}) = false
@inline endspace_correction(nspace::Int, count::Int) = count - nspace

@generated function Base.count{T<:Site,A<:Alphabet}(site::T, alg::BC, a::BioSequence{A}, b::BioSequence{A})
    n = bitsof(A)
    n_elems = div(64, n)

    quote
        if length(a) > length(b)
            return count(site, alg, b, a)
        end
        @assert length(a) ≤ length(b)

        nexta = bitindex(a, 1)
        nextb = bitindex(b, 1)
        stopa = bitindex(a, endof(a) + 1)
        counts = start_counter(site, a, b)
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
            counts = update_counter(counts, count_bitpar(T, A, x & m, y & m))
            #println("Count: ", counts)
            if correct_endspace(T, A)
                #println("N empty nibbles: ", nempty)
                nempty = $n_elems - div(k, $n)
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
                counts = update_counter(counts, count_bitpar(T, A, x, y))
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
                counts = update_counter(counts, count_bitpar(T, A, x & m, y & m))
                #println("Count: ", counts)
                if correct_endspace(T, A)
                    nempty = $n_elems - div(offs, $n)
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
                counts = update_counter(counts, count_bitpar(T, A, x, y))
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
                counts = update_counter(counts, count_bitpar(T, A, x & m, y & m))
                #println("Count: ", counts)
                if correct_endspace(T, A)
                    nempty = $n_elems - div(offs, $n)
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
