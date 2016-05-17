# Miscellaneous Utilities
# =======================
#
# Miscellaneous utilities for interactive tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
Copy a whole NucleotideSequence into the user's clipboard.

Useful for pasting a query sequence in web services like BLAST.
"""
function Base.clipboard(seq::BioSequence, width::Integer=50)
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
Base.clipboard(x::Kmer) = clipboard(string(x))
