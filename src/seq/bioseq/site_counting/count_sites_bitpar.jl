BC = BitparCount

@inline correct_emptyspace{T<:Site,A<:Alphabet}(::Type{T}, ::Type{A}) = false
@inline emptyspace_correction(nempty::Int, count::Int) = count - nempty

@generated function Base.count{S<:Site,A<:NucleicAcidAlphabets}(::Type{S}, ::Type{BC}, a::BioSequence{A}, b::BioSequence{A})
    n = bitsof(A)
    n_elems = div(64, n)

    quote
        if length(a) > length(b)
            return count(S, BC, b, a)
        end
        @assert length(a) ≤ length(b)

        nexta = bitindex(a, 1)
        nextb = bitindex(b, 1)
        stopa = bitindex(a, endof(a) + 1)
        counts = start_counter(S, A, A)
        spaces = 0

        # align reading position of `a.data` so that `offset(nexta) == 0`
        if nexta < stopa && offset(nexta) != 0
            x = a.data[index(nexta)] >> offset(nexta)
            y = b.data[index(nextb)] >> offset(nextb)
            if offset(nextb) > offset(nexta)
                y |= b.data[index(nextb)+1] << (64 - offset(nextb))
            end
            k = 64 - offset(nexta)
            m = mask(k)
            counts = update_counter(counts, count_bitpar(S, A, x & m, y & m))
            if correct_emptyspace(S, A)
                nempty = $n_elems - div(k, $n)
                counts = emptyspace_correction(nempty, counts)
            end
            nexta += k
            nextb += k
        end
        @assert offset(nexta) == 0

        if offset(nextb) == 0  # data are aligned with each other
            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                counts = update_counter(counts, count_bitpar(S, A, x, y))
                nexta += 64
                nextb += 64
            end

            if nexta < stopa
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                offs = stopa - nexta
                m = mask(offs)
                counts = update_counter(counts, count_bitpar(S, A, x & m, y & m))
                if correct_emptyspace(S, A)
                    nempty = $n_elems - div(offs, $n)
                    counts = emptyspace_correction(nempty, counts)
                end
            end
        elseif nexta < stopa
            y = b.data[index(nextb)]
            nextb += 64

            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                z = b.data[index(nextb)]
                y = y >> offset(nextb) | z << (64 - offset(nextb))

                counts = update_counter(counts, count_bitpar(S, A, x, y))
                y = z
                nexta += 64
                nextb += 64
            end

            if nexta < stopa
                x = a.data[index(nexta)]
                y = y >> offset(nextb)
                if 64 - offset(nextb) < stopa - nexta
                    y |= a.data[index(nextb)] << (64 - offset(nextb))
                end
                offs = stopa - nexta
                m = mask(offs)
                counts = update_counter(counts, count_bitpar(S, A, x & m, y & m))
                if correct_emptyspace(S, A)
                    nempty = $n_elems - div(offs, $n)
                    counts = emptyspace_correction(nempty, counts)
                end
            end
        end
        return counts
    end
end
