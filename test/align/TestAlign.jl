module TestAlign

using FactCheck
using Bio.Align

facts("Alignments") do
    context("Operations") do
        context("Conversions") do
            ops = [OP_GAP, OP_MATCH, OP_MM, OP_MM, OP_N, OP_N,
            OP_MISMATCH, OP_MISMATCH, OP_SCLIP, OP_SCLIP,
            OP_HCLIP, OP_HCLIP, OP_INSERT, OP_INSERT,
            OP_DELETE, OP_DELETE, OP_PAD, OP_PAD]
            @inbounds for i in 1:length(CIGAR_CHARS)
                @fact Operation(CIGAR_CHARS[i]) --> ops[i]
            end
            ops = unique(ops)
            chars = ['-', '=', 'M', 'N', 'X', 'S', 'H', 'I', 'D', 'P']
            @inbounds for i in 1:length(ops)
                @fact Char(eval(ops[i])) --> chars[i]
            end
        end
    end

    context("CIGAR") do

    end

    context("Anchors") do

    end
end

end # TestAlign
