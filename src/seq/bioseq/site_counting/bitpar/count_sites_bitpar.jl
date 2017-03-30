
"""
    count_sites_bitpar()
"""
@generated function count_sites_bitpar{T<:Site,A}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})
    n = bitsof(A)
    if n == 2
        count_func = :(bitpar_mismatches2)
    elseif n == 4
        count_func = :count_nibbles
    else
        error("n (= $n) ∉ (2, 4)")
    end

    quote
        if length(a) > length(b)
            return wee(b, a)
        end
        @assert length(a) ≤ length(b)

        nexta = bitindex(a, 1)
        nextb = bitindex(b, 1)
        stopa = bitindex(a, endof(a) + 1)
        counts = 0

        # align reading position of `a.data` so that `offset(nexta) == 0`
        if nexta < stopa && offset(nexta) != 0
            println("Aligning sequence A so it has no offset.")
            x = a.data[index(nexta)] >> offset(nexta)
            y = b.data[index(nextb)] >> offset(nextb)
            if offset(nextb) > offset(nexta)
                y |= b.data[index(nextb)+1] << (64 - offset(nextb))
            end
            k = 64 - offset(nexta)
            m = mask(k)
            println("mask: $(hex(m)), k: $(hex(k))")
            println("x: $(hex(x)), y: $(hex(y))")
            println("x masked: $(hex(x & m)), y masked: $(hex(y & m))")
            println("x type: $(typeof(x))")
            counts += $count_func(T, x & m, y & m)
            nexta += k
            nextb += k
            println("Done initial alignment.")
        end
        @assert offset(nexta) == 0

        if offset(nextb) == 0  # data are aligned with each other
            println("A and B are aligned with each other.")
            println("Doing aligned loop.")
            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                println("x: $(hex(x))")
                println("y: $(hex(y))")
                counts += $count_func(T, x, y)
                nexta += 64
                nextb += 64
            end
            println("Done aligned loop.")

            if nexta < stopa
                println("Doing aligned tail.")
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                m = mask(stopa - nexta)
                println("x: $(hex(x))")
                println("y: $(hex(y))")
                println("m: $(hex(m))")
                println("x masked: $(hex(x & m))")
                println("y masked: $(hex(y & m))")
                counts += $count_func(T, x & m, y & m)
                println("Done aligned tail.")
            end
        elseif nexta < stopa
            println("A and B are not aligned.")
            println("Doing unaligned loop.")

            y = b.data[index(nextb)]
            nextb += 64

            println("stopa: $stopa, nexta: $nexta")

            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                z = b.data[index(nextb)]
                y = y >> offset(nextb) | z << (64 - offset(nextb))

                println("x: $(hex(x)), z: $(hex(z)), y: $(hex(y))")
                println("x type: $(typeof(x))")
                counts += $count_func(T, x, y)
                y = z
                nexta += 64
                nextb += 64
            end
            println("Done unaligned loop.")

            if nexta < stopa
                println("Doing unaligned tail.")
                x = a.data[index(nexta)]
                y = y >> offset(nextb)
                println("x: $(hex(x)), y: $(hex(y))")
                println("x type: $(typeof(x))")
                if 64 - offset(nextb) < stopa - nexta
                    y |= a.data[index(nextb)] << (64 - offset(nextb))
                end
                println("y: $(hex(y))")
                m = mask(stopa - nexta)
                println("x masked: $(hex(x & m)), y masked: $(hex(y & m))")
                println("x type: $(typeof(x))")
                counts += $count_func(T, x & m, y & m)
                println("Done unaligned tail.")
            end
        end

        println("All done.")
        return counts
    end
end
