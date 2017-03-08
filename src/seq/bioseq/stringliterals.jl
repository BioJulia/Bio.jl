# String Decorators
# -----------------

remove_newlines(s) = replace(s, r"\r|\n", "")

macro dna_str(seq, flag)
    if flag == "s"
        return DNASequence(remove_newlines(seq))
    elseif flag == "d"
        return quote
            DNASequence($(remove_newlines(seq)))
        end
    end
    error("Invalid DNA flag: '$(flag)'")
end

macro dna_str(seq)
    return DNASequence(remove_newlines(seq))
end

macro rna_str(seq, flag)
    if flag == "s"
        return RNASequence(remove_newlines(seq))
    elseif flag == "d"
        return quote
            RNASequence($(remove_newlines(seq)))
        end
    end
    error("Invalid RNA flag: '$(flag)'")
end

macro rna_str(seq)
    return RNASequence(remove_newlines(seq))
end

macro aa_str(seq, flag)
    if flag == "s"
        return AminoAcidSequence(remove_newlines(seq))
    elseif flag == "d"
        return quote
            AminoAcidSequence($(remove_newlines(seq)))
        end
    end
    error("Invalid Amino Acid flag: '$(flag)'")
end

macro aa_str(seq)
    return AminoAcidSequence(remove_newlines(seq))
end

macro char_str(seq, flag)
    if flag == "s"
        return CharSequence(remove_newlines(seq))
    elseif flag == "d"
        return quote
            CharSequence($(remove_newlines(seq)))
        end
    end
    error("Invalid Char flag: '$(flag)'")
end

macro char_str(seq)
    return CharSequence(remove_newlines(seq))
end
