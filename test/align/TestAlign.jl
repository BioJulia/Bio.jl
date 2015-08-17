module TestAlign

using FactCheck
using Bio.Align

chars1 = ['-', '=', 'M', 'm', 'N', 'n', 'X', 'x', 'S', 's', 'H', 'h', 'I', 'i',
         'D', 'd', 'P', 'p']

chars2 = ['-', '=', 'M', 'N', 'X', 'S', 'H', 'I', 'D', 'P']

ops1 = [OP_GAP, OP_MATCH, OP_MM, OP_MM, OP_N, OP_N, OP_MISMATCH, OP_MISMATCH,
       OP_SCLIP, OP_SCLIP, OP_HCLIP, OP_HCLIP, OP_INSERT, OP_INSERT, OP_DELETE,
       OP_DELETE, OP_PAD, OP_PAD]

ops2 = unique(ops1)

cigarSize = abs(rand(Int, length(ops2)))

strings = ["$(cigarSize[i])$(chars2[i])" for i in 1:length(cigarSize)]

function cigarsMatch(cigar::CIGAR, op::Operation, size::Int)
    cigar.OP == op && cigar.Size == size
end

function randomCIGARS()

end


facts("Alignments") do
    context("Operations") do
        context("Conversions") do
            @inbounds for i in 1:length(chars1)
                @fact Operation(chars1[i]) --> ops1[i]
            end
            @inbounds for i in 1:length(ops2)
                @fact Char(eval(ops2[i])) --> chars2[i]
            end
        end

    end

    context("CIGAR") do
        context("Constructors") do
            context("From Operation and Size") do
                @inbounds for i in 1:length(ops2)
                    @fact cigarsMatch(CIGAR(ops2[i], cigarSize[i]), ops2[i], cigarSize[i]) --> true
                end
            end

            context("From Character and Size") do
                @inbounds for i in 1:length(chars2)
                    @fact cigarsMatch(CIGAR(chars2[i], cigarSize[i]), ops2[i], cigarSize[i]) --> true
                end
            end

            context("From Strings") do
                @inbounds for i in 1:length(strings)
                    @fact cigarsMatch(CIGAR(strings[i])) --> true
            end





        end
    end

    context("Anchors") do

    end
end

end # TestAlign
