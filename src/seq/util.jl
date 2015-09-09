"""
Copy a whole NucleotideSequence into the user's clipboard.

Useful for pasting a query sequence in web services like BLAST.
"""
function clipboard(seq::Union(NucleotideSequence,AminoAcidSequence), width::Integer=50)
    buf = IOBuffer()
    for i in 1:length(seq)
        if i % width == 1 && i > width
            write(buf, '\n')
        end
        write(buf, convert(Char, seq[i]))
    end
    clipboard(bytestring(buf))
end

"""
Copy a Kmer into the user's clipboard.

Useful for pasting a query sequence in web services like BLAST.
"""
clipboard(x::Kmer) = clipboard(convert(AbstractString, x))
