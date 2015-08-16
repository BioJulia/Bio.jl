module TestAlign

using FactCheck
using Bio.Align

facts("Alignments") do
    context("Operations") do
        context("Conversions") do
            context("From characters to Operations") do
                chars = "-=MmNnXxSsHhIiDdPp"
                ops = [OP_GAP, OP_MATCH, OP_MM, OP_MM, OP_N, OP_N,
                OP_MISMATCH, OP_MISMATCH, OP_SCLIP, OP_SCLIP,
                OP_HCLIP, OP_HCLIP, OP_INSERT, OP_INSERT,
                OP_DELETE, OP_DELETE, OP_PAD, OP_PAD]
                @inbounds for i in 1:length(chars)
                    @fact Operation(chars[i]) --> ops[i] "Conversion from $(chars[i]) to $(ops[i])"
                end
            end
        end
    end

    context("CIGAR") do

    end

    context("Anchors") do

    end
end

end # TestAlign
