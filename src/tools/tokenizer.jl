module Tokenizer

using Base.Intrinsics

# ========================
# A Tokenizer for Strings
# ========================


@Docile.doc """
Tokenizer type that is responsible for splitting a string into token according to a
regex specification. The regex specification is stored in `dict`.

**Fields:**

* `dict`:      A Dictionary of String keys and Regex values.
* `tokenizer`: A regex that is generated from the specification.
""" ->
type Tokenizer
  dict::Dict{String, Regex}
  tokenizer::Regex
  """
  Inner constructor for Tokenizers. Accepts a vector of String, Regex
  tuples that makes the token specifications.

  **Parameters:**

  * `x`: An array of ASCIIString, Regex tuples.

  **Returns:** Instance of `Tokenizer`
  """
  function Tokenizer(x::Vector{Tuple{ASCIIString, Regex}})
    dictvalues = [i[2] for i in x]
    stringvalues = String[i.pattern for i in dictvalues]
    combinedstring = join(stringvalues, "|")
    finalstring = "($combinedstring)"
    return new([i => j for (i, j) in x], Regex(finalstring))
  end
end


end
