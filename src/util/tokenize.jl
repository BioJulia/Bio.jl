# util/tokenize.jl
# ================
#
# A module defining types and methods for tokenizing strings.
#
# Useful in creating parsers.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md


module Tokenize

using Compat
import Compat.String

export
    Tokenizer, TokenizerResult, tokenize


"""
Tokenizer type that is responsible for splitting a string into token according to a
regex specification. The regex specification is stored in `dict`.

**Fields:**

* `dict`:      A Dictionary of String keys and Regex values.
* `tokenizer`: A regex that is generated from the specification.
"""
type Tokenizer
    dict::Dict{String, Regex}
    tokenizer::Regex
    """
    Inner constructor for Tokenizers. Accepts a vector of String, Regex
    tuples that makes the token specifications.

    **Parameters:**

    * `x`: An array of `Compat.String`, `Regex` tuples.

    **Returns:** Instance of `Tokenizer`
    """
    function Tokenizer(x::Vector{Tuple{String, Regex}})
        dictvalues = [i[2] for i in x]
        stringvalues = Compat.String[i.pattern for i in dictvalues]
        combinedstring = join(stringvalues, "|")
        finalstring = "($combinedstring)"
        return new([i => j for (i, j) in x], Regex(finalstring))
    end
end

immutable TokenizerResult
    tokens::Vector{String}
    failiures::Vector{String}

    function TokenizerResult(tokens = String[], failiures = String[])
        return new(tokens, failiures)
    end
end

function tokenize(s::String, t::Tokenizer)
    lastpos = 0
    out = TokenizerResult()
    for mat in eachmatch(t.tokenizer, s)
        if mat.offset > lastpos + 1
            push!(failiures, str[lastpos + 1 : mat.offset - 1])
        end
        push!(out.tokens, mat.match)
        lastpos = mat.offset + length(mat.match) - 1
    end
    return out
end

end
