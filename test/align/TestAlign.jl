module TestAlign

using FactCheck
using Bio.Align

## Constants and Parameters used for defining and running the tests.

const MAX_CIGAR_SIZE = 100000

const NUMBER_OF_RANDOM_CIGAR = 10

const MAX_CIGARS_IN_CIGAR_STRING = 50

const NUMBER_OF_RANDOM_CIGAR_STRINGS = 10

const ALLOWED_CHARS = ['-', '=', 'M', 'm', 'N', 'n', 'X', 'x', 'S', 's', 'H',
                       'h', 'I', 'i', 'D', 'd', 'P', 'p']

const CHARS_FROM_OPS = ['-', '=', 'M', 'N', 'X', 'S', 'H', 'I', 'D', 'P']

const OPS_TWICE = [OP_GAP, OP_MATCH, OP_MM, OP_MM, OP_N, OP_N, OP_MISMATCH,
              OP_MISMATCH, OP_SCLIP, OP_SCLIP, OP_HCLIP, OP_HCLIP, OP_INSERT,
              OP_INSERT, OP_DELETE, OP_DELETE, OP_PAD, OP_PAD]

const OPS_UNIQUE = unique(OPS_TWICE)

# Need some random CIGAR sizes for some tests.
const cigarSize = abs(rand(1:MAX_CIGAR_SIZE, length(OPS_UNIQUE)))

strings = ["$(cigarSize[i])$(CHARS_FROM_OPS[i])" for i in 1:length(cigarSize)]

## Functions used for tests.

function cigarsMatch(cigar::CIGAR, op::Operation, size::Int)
    cigar.OP == op && cigar.Size == size
end


## Test procedure.

facts("Alignments") do
    context("Operations") do
        context("Constructors and Conversions") do
            @inbounds for i in 1:length(chars1)
                @fact Operation(chars1[i]) --> OPS_TWICE[i]
            end
            @inbounds for i in 1:length(OPS_UNIQUE)
                @fact Char(eval(OPS_UNIQUE[i])) --> CHARS_FROM_OPS[i]
            end
        end

    end

    context("Single CIGARs") do
        context("Constructors and Conversions") do
            context("From Operation and Size") do
                # Check GICARs are correctly constructed from an Operation and Size argument.
                @inbounds for i in 1:length(OPS_UNIQUE)
                    @fact cigarsMatch(CIGAR(OPS_UNIQUE[i], cigarSize[i]), OPS_UNIQUE[i], cigarSize[i]) --> true
                end
            end

            context("From Character and Size") do
                # Check CIGARs are correctly constructed from a character and size argument.
                @inbounds for i in 1:length(CHARS_FROM_OPS)
                    @fact cigarsMatch(CIGAR(CHARS_FROM_OPS[i], cigarSize[i]), OPS_UNIQUE[i], cigarSize[i]) --> true
                end
            end

            context("From Strings") do
                # Check Strings are successfully converted to CIGAR.
                @inbounds for i in 1:length(strings)
                    @fact cigarsMatch(CIGAR(strings[i]), OPS_UNIQUE[i], cigarSize[i]) --> true
                end
            end

            context("To Strings from CIGAR") do
                # Check that CIGAR operations are successfully converted to Strings.
                @inbounds for i in 1:length(OPS_UNIQUE)
                    @fact String(CIGAR(OPS_UNIQUE[i], cigarSize[i])) --> "$(cigarSize[i])$(OPS_UNIQUE[i])"
                end
            end
        end
    end

    context("CIGARS") do
        context("Constructors and Conversions") do
            context("From String to CIGARS") do
                for i in 1:NUMBER_OF_RANDOM_CIGAR_STRINGS
                    numCigar = rand(1:MAX_CIGARS_IN_CIGAR_STRING)
                    operationChars = rand(ALLOWED_CHARS, numCigar)
                    operationOPS = [Operation(i) for i in operationChars]
                    cigarLengths = rand(1:MAX_CIGAR_SIZE, numCigar)
                    out = ""
                    @inbounds for i in 1:MAX_CIGARS_IN_CIGAR_STRING
                        out *= "$(cigarLengths[i])$(operationChars[i])"
                    end
                    cigars = CIGARS(out)
                    @fact length(cigars) --> numCigar
                    for i in 1:numCigar
                        @fact cigars[i].OP --> operationOPS[i]
                        @fact cigars[i].Size --> cigarLengths[i]
                    end
                end
            end
        end
    end

    context("AlignmentAnchor") do

    end
end

end # TestAlign
