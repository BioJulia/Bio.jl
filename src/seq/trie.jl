type SeqTrie{A<:Alphabet, T}
    value::T
    children::Dict{A,SeqTrie{A, T}}
    is_key::Bool

    function SeqTrie()
        self = new()
        self.children = Dict{A,SeqTrie{A, T}}()
        self.is_key = false
        return self
    end

    function SeqTrie(kv)
        t = SeqTrie{A, T}()
        for (k,v) in kv
            t[k] = v
        end
        return t
    end
    SeqTrie(ks, vs) = SeqTrie(zip(ks, vs))
end

# Extra ctors
## Typed values, aliased keys
DNATrie{T}() = SeqTrie{DNAAlphabet, T}()
RNATrie{T}() = SeqTrie{RNAAlphabet, T}()
AminoAcidTrie{T}() = SeqTrie{AminoAcidAlphabet, T}()

## Any values, aliased keys
SeqTrie{A}() = SeqTrie{A, Any}()
DNATrie() = SeqTrie{DNAAlphabet, Any}()
RNATrie() = SeqTrie{RNAAlphabet, Any}()
AminoAcidTrie() = SeqTrie{AminoAcidAlphabet, Any}()

function setindex!{A, T}(t::SeqTrie{A, T}, val::T, key::BioSequence{A})
    node = t
    for char in key
        if !haskey(node.children, char)
            node.children[char] = SeqTrie{T}()
        end
        node = node.children[char]
    end
    node.is_key = true
    node.value = val
end

function getindex{A, T}(t::SeqTrie{A, T}, key::BioSequence{A})
    node = subtrie(t, key)
    if node != nothing && node.is_key
        return node.value
    end
    throw(KeyError("key not found: $key"))
end

function subtrie{A, T}(t::SeqTrie{A, T}, prefix::BioSequence{A})
    node = t
    for char in prefix
        if !haskey(node.children, char)
            return nothing
        else
            node = node.children[char]
        end
    end
    node
end

function haskey{A, T}(t::SeqTrie{A, T}, key::BioSequence{A})
    node = subtrie(t, key)
    return node != nothing && node.is_key
end

function hasprefix{A, T}(t::SeqTrie{A, T}, prefix::BioSequence{A})
    node = subtrie(t, key)
    return node != nothing
end

function get{A, T}(t::SeqTrie{A, T}, key::BioSequence{A}, notfound)
    node = subtrie(t, key)
    if node != nothing && node.is_key
        return node.value
    end
    notfound
end

function keys{A, T}(t::SeqTrie{A, T}, prefix::BioSequence{A}=BioSequence{A}(),
                    found=BioSequence{A}[])
    if t.is_key
        push!(found, prefix)
    end
    for (char,child) in t.children
        newprefix = prefix * BioSequence{A}(char)
        keys(child, newprefix, found)
    end
    return found
end

function keys_with_prefix{A, T}(t::SeqTrie{A, T}, prefix::BioSequence{A})
    st = subtrie(t, prefix)
    st != nothing ? keys(st,prefix) : []
end

# The state of a SeqTrieIterator is a pair (t::SeqTrie, i::Int),
# where t is the SeqTrie which was the output of the previous iteration
# and i is the index of the current character of the string.
# The indexing is potentially confusing;
# see the comments and implementation below for details.
immutable SeqTrieIterator{A, T}
    t::SeqTrie{A, T}
    seq::BioSequence{A}
end

# At the start, there is no previous iteration,
# so the first element of the state is undefined.
# We use a "dummy value" of it.t to keep the type of the state stable.
# The second element is 0
# since the root of the trie corresponds to a length 0 prefix of seq.
start(it::SeqTrieIterator) = (it.t, 0)

function next(it::SeqTrieIterator, state)
    t, i = state
    i == 0 && return it.t, (it.t, 1)

    t = t.children[it.seq[i]]
    return (t, (t, i + 1))
end

function done(it::SeqTrieIterator, state)
    t, i = state
    i == 0 && return false
    i == length(it.seq) + 1 && return true
    return !(it.seq[i] in keys(t.children))
end

path(t::SeqTrie, seq::BioSequence) = SeqTrieIterator(t, seq)
if VERSION >= v"0.5.0-dev+3294"
    Base.iteratorsize(::Type{SeqTrieIterator}) = Base.SizeUnknown()
end
