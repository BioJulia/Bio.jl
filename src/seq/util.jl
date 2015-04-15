# Copy a whole sequence into the user's clipboard; it would be useful for
# pasting a query sequence in web services like BLAST
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

clipboard(x::Kmer) = clipboard(convert(String, x))
