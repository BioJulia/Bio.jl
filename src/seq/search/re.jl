# Regular Expression
# ==================
#
# Regular expression sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# String Decorators
# -----------------

macro biore_str(pat, opt...)
    if isempty(opt)
        error("symbol option is required: d(na), r(na), or a(a)")
    end
    @assert length(opt) == 1
    opt = opt[1]

    if opt ∈ ("d", "dna")
        :($(RE.Regex{DNANucleotide}(pat, :pcre)))
    elseif opt ∈ ("r", "rna")
        :($(RE.Regex{RNANucleotide}(pat, :pcre)))
    elseif opt ∈ ("a", "aa")
        :($(RE.Regex{AminoAcid}(pat, :pcre)))
    else
        error("invalid symbol option: $(opt)")
    end
end

macro prosite_str(pat)
    :($(RE.Regex{AminoAcid}(pat, :prosite)))
end


module RE

using ..Seq

# Syntax tree
# -----------

type SyntaxTree
    head::Symbol
    args::Vector{Any}

    function SyntaxTree(head, args)
        @assert head ∈ (
            :|,
            :*, :+, Symbol("?"), :range,
            Symbol("*?"), Symbol("+?"), Symbol("??"), Symbol("range?"),
            :set, :compset, :sym, :bits,
            :capture, :nocapture, :concat, :head, :last)
        return new(head, args)
    end
end

expr(head, args) = SyntaxTree(head, args)

function charset(s)
    set = Set{Char}()
    union!(set, s)
    union!(set, lowercase(s))
    return set
end

# list of symbols available for each symbol type
const symbols = ObjectIdDict(
    DNANucleotide => charset("ACGTMRWSYKVHDBN"),
    RNANucleotide => charset("ACGUMRWSYKVHDBN"),
    AminoAcid     => charset("ARNDCQEGHILKMFPSTWYVOUBJZX"))

macro check(ex, err)
    quote
        if !$ex
            throw($err)
        end
    end
end

function parse{T}(::Type{T}, pat::AbstractString)
    parens = Char[]  # stack of parens
    ex, _ = parserec(T, pat, start(pat), parens)
    @check isempty(parens) ArgumentError("'(' is not closed")
    return ex
end

function parserec{T}(::Type{T}, pat, s, parens)
    args = []
    while !done(pat, s)
        c, s = next(pat, s)
        if c == '*'
            @check !isempty(args) ArgumentError("unexpected '*'")
            arg = pop!(args)
            push!(args, expr(:*, [arg]))
        elseif c == '+'
            @check !isempty(args) ArgumentError("unexpected '+'")
            arg = pop!(args)
            push!(args, expr(:+, [arg]))
        elseif c == '?'
            @check !isempty(args) ArgumentError("unexpected '?'")
            arg = pop!(args)
            if arg.head ∈ (:*, :+, Symbol("?"), :range)
                # lazy quantifier
                push!(args, expr(Symbol(arg.head, '?'), arg.args))
            else
                # zero-or-one quantifier
                push!(args, expr(Symbol("?"), [arg]))
            end
        elseif c == '{'
            @check !isempty(args) ArgumentError("unexpected '{'")
            rng, s = parserange(pat , s)
            arg = pop!(args)
            push!(args, expr(:range, [rng, arg]))
        elseif c == '|'
            arg1 = expr(:concat, args)
            arg2, s = parserec(T, pat, s, parens)
            args = []
            push!(args, expr(:|, [arg1, arg2]))
        elseif c == '['
            setexpr, s = parseset(T, pat, s)
            push!(args, setexpr)
        elseif c == '('
            push!(parens, '(')
            head = :capture
            if peek(pat, s) == '?'
                _, s = next(pat, s)
                if peek(pat, s) == ':'
                    head = :nocapture
                    _, s = next(pat, s)
                end
            end
            arg, s = parserec(T, pat, s, parens)
            push!(args, expr(head, [arg]))
        elseif c == ')'
            @check !isempty(parens) && parens[end] == '(' ArgumentError("unexpected ')'")
            pop!(parens)
            return expr(:concat, args), s
        elseif c == '^'
            push!(args, expr(:head, []))
        elseif c == '$'
            push!(args, expr(:last, []))
        elseif isspace(c)
            # skip
        elseif c ∈ symbols[T]
            push!(args, expr(:sym, [convert(T, c)]))
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    return expr(:concat, args), s
end

function parserange(pat, s)
    lo = hi = -1
    comma = false
    while !done(pat, s)
        c, s = next(pat, s)
        if isdigit(c)
            d = c - '0'
            if comma
                if hi < 0
                    hi = 0
                end
                hi = 10hi + d
            else
                if lo < 0
                    lo = 0
                end
                lo = 10lo + d
            end
        elseif c == ','
            comma = true
        elseif c == '}'
            break
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    if comma
        if lo < 0 && hi < 0
            throw(ArgumentError("invalid range"))
        elseif lo ≥ 0 && hi < 0
            return (lo,), s
        elseif lo < 0 && hi ≥ 0
            return (0, hi), s
        else  # lo ≥ 0 && hi ≥ 0
            if lo > hi
                throw(ArgumentError("invalid range"))
            end
            return (lo, hi), s
        end
    else
        return lo, s
    end
end

function peek(pat, s)
    if done(pat, s)
        throw(ArgumentError("unexpected end of pattern"))
    end
    return next(pat, s)[1]
end

function parseset{T}(::Type{T}, pat, s)
    if peek(pat, s) == '^'
        head = :compset
        _, s = next(pat, s)
    else
        head = :set
    end
    set = T[]
    while !done(pat, s)
        c, s = next(pat, s)
        if c ∈ symbols[T]
            push!(set, convert(T, c))
        elseif c == ']'
            break
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    return expr(head, set), s
end

function parse_prosite(pat)
    s = start(pat)
    args = []
    while !done(pat, s)
        c, s = next(pat, s)
        if c == '['
            set, s = parseset_prosite(pat, s, ']')
            push!(args, set)
        elseif c == '{'
            set, s = parseset_prosite(pat, s, '}')
            push!(args, set)
        elseif c == '('
            @check !isempty(args) ArgumentError("unexpected '('")
            rng, s = parserange_prosite(pat, s)
            arg = pop!(args)
            push!(args, expr(:range, [rng, arg]))
        elseif c == '-'
            # concat
            continue
        elseif c == 'x'
            push!(args, expr(:sym, [AA_X]))
        elseif c == '<'
            push!(args, expr(:head, []))
        elseif c == '>'
            push!(args, expr(:last, []))
        elseif c ∈ symbols[AminoAcid]
            push!(args, expr(:sym, [convert(AminoAcid, c)]))
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    return expr(:concat, args)
end

function parserange_prosite(pat, s)
    lo = hi = -1
    comma = false
    while !done(pat, s)
        c, s = next(pat, s)
        if isdigit(c)
            d = c - '0'
            if comma
                if hi < 0
                    hi = 0
                end
                hi = 10hi + d
            else
                if lo < 0
                    lo = 0
                end
                lo = 10lo + d
            end
        elseif c == ','
            comma = true
        elseif c == ')'
            break
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    if comma
        if lo < 0 || hi < 0
            throw(ArgumentError("invalid range"))
        end
        return (lo, hi), s
    else
        if lo < 0
            throw(ArgumentError("invalid range"))
        end
        return lo, s
    end
end

function parseset_prosite(pat, s, close)
    set = AminoAcid[]
    while !done(pat, s)
        c, s = next(pat, s)
        if c ∈ symbols[AminoAcid]
            push!(set, convert(AminoAcid, c))
        elseif c == close
            if close == ']'
                return expr(:set, set), s
            elseif close == '}'
                return expr(:compset, set), s
            end
            @assert false
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
end

function bits2sym{T}(::Type{T}, bits::UInt32)
    for x in alphabet(T)
        if Seq.compatbits(x) == bits
            return x
        end
    end
    error("bits are not found")
end

mask{T<:Nucleotide}(::Type{T}) = (UInt32(1) << 4) - one(UInt32)
@assert Int(AA_U) + 1 == 22  # check there are 22 unambiguous amino acids
mask(::Type{AminoAcid}) = (UInt32(1) << 22) - one(UInt32)

function desugar{T}(::Type{T}, tree::SyntaxTree)
    head = tree.head
    args = tree.args
    if head == :+
        # e+ => ee*
        head = :concat
        args = [args[1], expr(:*, [args[1]])]
    elseif head == Symbol("+?")
        # e+? => ee*?
        head = :concat
        args = [args[1], expr(Symbol("*?"), [args[1]])]
    elseif head == Symbol("?")
        # e? => e|
        head = :|
        args = [args[1], expr(:concat, [])]
    elseif head == Symbol("??")
        # e?? => |e
        head = :|
        args = [expr(:concat, []), args[1]]
    elseif head == :sym
        head = :bits
        args = [Seq.compatbits(args[1])]
    elseif head == :set
        bits = UInt32(0)
        for arg in args
            bits |= Seq.compatbits(arg)
        end
        head = :bits
        args = [bits]
    elseif head == :compset
        bits = UInt32(0)
        for arg in args
            bits |= Seq.compatbits(arg)
        end
        head = :bits
        args = [~bits & mask(T)]
    elseif head == :nocapture
        head = :concat
    elseif head == :range || head == Symbol("range?")
        rng = args[1]
        pat = args[2]
        greedy = head == :range
        if isa(rng, Int)
            # e{m} => eee...e
            #         |<-m->|
            head = :concat
            args = replicate(pat, rng)
        elseif isa(rng, Tuple{Int})
            # e{m,} => eee...ee*   (greedy)
            #          |<-m->|
            # e{m,} => eee...ee*?  (lazy)
            #          |<-m->|
            head = :concat
            args = replicate(pat, rng[1])
            if greedy
                push!(args, expr(:*, [pat]))
            else
                push!(args, expr(Symbol("*?"), [pat]))
            end
        elseif isa(rng, Tuple{Int,Int})
            # e{m,n} => eee...eee|eee...ee|...|eee...e  (greedy)
            #           |<- n ->|              |<-m->|
            # e{m,n} => eee...e|eee...ee|...|eee...eee  (lazy)
            #           |<-m->|              |<- n ->|
            m, n = rng
            @assert m ≤ n
            if m == n
                head = :concat
                args = replicate(pat, m)
            else
                head = :|
                f = (k, t) -> expr(:|, [expr(:concat, replicate(pat, k)), t])
                if greedy
                    args = foldr(f, expr(:concat, replicate(pat, m)), n:-1:m+1).args
                else
                    args = foldr(f, expr(:concat, replicate(pat, n)), m:1:n-1).args
                end
            end
        else
            @assert false "invalid AST"
        end
    end
    i = 1
    for i in 1:endof(args)
        args[i] = desugar(T, args[i])
    end
    return expr(head, args)
end

function desugar{T}(::Type{T}, atom)
    return atom
end

function replicate(x, n)
    return collect(repeated(x, n))
end


# Compiler
# --------

#  tag  | operation |           meaning
# ------|-----------|----------------------------
# 0b000 | match     | the pattern matches
# 0b001 | bits b    | matches b in bit-wise way
# 0b010 | jump l    | jump to l
# 0b011 | push l    | push l and go to next
# 0b100 | save i    | save string state to i
# 0b101 | head      | matches the head of string
# 0b110 | last      | matches the last of string
# 0b111 | fork l    | push next and go to l

bitstype 32 Op

const MatchTag = UInt32(0b000) << 29
const BitsTag  = UInt32(0b001) << 29
const JumpTag  = UInt32(0b010) << 29
const PushTag  = UInt32(0b011) << 29
const SaveTag  = UInt32(0b100) << 29
const HeadTag  = UInt32(0b101) << 29
const LastTag  = UInt32(0b110) << 29
const ForkTag  = UInt32(0b111) << 29

# constructors
match() = reinterpret(Op, MatchTag)
bits(b::UInt32) = reinterpret(Op, BitsTag | b)
jump(l::Int) = reinterpret(Op, JumpTag | UInt32(l))
push(l::Int) = reinterpret(Op, PushTag | UInt32(l))
save(l::Int) = reinterpret(Op, SaveTag | UInt32(l))
head() = reinterpret(Op, HeadTag)
last() = reinterpret(Op, LastTag)
fork(l::Int) = reinterpret(Op, ForkTag | UInt32(l))

const operand_mask = (UInt32(1) << 29) - one(UInt32)

function tag(op::Op)
    return reinterpret(UInt32, op) & ~operand_mask
end

function operand(op::Op)
    return reinterpret(UInt32, op) &  operand_mask
end

function Base.show(io::IO, op::Op)
    t = tag(op)
    x = operand(op)
    if t == MatchTag
        print(io, "match")
    elseif t == BitsTag
        print(io, "bits ", @sprintf("0x%08x", x))
    elseif t == JumpTag
        print(io, "jump ", x)
    elseif t == PushTag
        print(io, "push ", x)
    elseif t == SaveTag
        print(io, "save ", x)
    elseif t == HeadTag
        print(io, "head")
    elseif t == LastTag
        print(io, "last")
    elseif t == ForkTag
        print(io, "fork ", x)
    else
        @assert false
    end
end

function print_program(prog::Vector{Op})
    L = endof(prog)
    for l in 1:L
        print(lpad(l, ndigits(L)), ": ", prog[l])
        if l != endof(prog)
            println()
        end
    end
end

function compile(tree::SyntaxTree)
    code = Op[]
    push!(code, save(1))
    compilerec!(code, tree, 2)
    push!(code, save(2))
    push!(code, match())
    return code
end

function compilerec!(code, tree::SyntaxTree, k)
    h = tree.head
    args = tree.args
    if h == :bits
        push!(code, bits(UInt32(args[1])))
    elseif h == :concat
        for arg in args
            k = compilerec!(code, arg, k)
        end
    elseif h == :*
        push!(code, push(0))  # placeholder
        l = length(code)
        k = compilerec!(code, args[1], k)
        push!(code, jump(l))
        code[l] = push(length(code) + 1)
    elseif h == Symbol("*?")
        push!(code, fork(0))  # placeholder
        l = length(code)
        k = compilerec!(code, args[1], k)
        push!(code, jump(l))
        code[l] = fork(length(code) + 1)
    elseif h == :|
        push!(code, push(0))  # placeholder
        l = length(code)
        k = compilerec!(code, args[1], k)
        push!(code, jump(0))  # placeholder
        code[l] = push(length(code) + 1)
        l = length(code)
        k = compilerec!(code, args[2], k)
        code[l] = jump(length(code) + 1)
    elseif h == :capture
        k′ = k
        push!(code, save(2k′-1))
        k = compilerec!(code, args[1], k + 1)
        push!(code, save(2k′))
    elseif h == :head
        push!(code, head())
    elseif h == :last
        push!(code, last())
    else
        @assert false "invalid tree"
    end
    return k
end


# Virtual machine
# ---------------

"""
Regular expression for `DNANucleotide`, `RNANucleotide`, and `AminoAcid`.
"""
immutable Regex{T}
    pat::String       # regular expression pattern (for printing)
    code::Vector{Op}  # compiled code
    nsaves::Int       # the number of `save` operations in `code`

    function Regex(pat::AbstractString, syntax=:pcre)
        if syntax == :pcre
            ast = desugar(T, parse(T, pat))
        elseif syntax == :prosite
            if T != AminoAcid
                throw(ArgumentError("alphabet must be AminoAcid for PROSITE syntax"))
            end
            ast = desugar(AminoAcid, parse_prosite(pat))
        else
            throw(ArgumentError("invalid syntax: $syntax"))
        end
        code = compile(ast)
        nsaves = 0
        for op in code
            nsaves += tag(op) == SaveTag
        end
        @assert iseven(nsaves)
        return new(pat, code, nsaves)
    end
end

function Base.show{T}(io::IO, re::Regex{T})
    if T == DNANucleotide
        opt = "dna"
    elseif T == RNANucleotide
        opt = "rna"
    elseif T == AminoAcid
        opt = "aa"
    else
        assert(false)
    end
    print(io, "biore\"", re.pat, "\"", opt)
end

"""
Result of matching by `Regex`.
"""
immutable RegexMatch{S}
    seq::S
    captured::Vector{Int}
end

function Base.show(io::IO, m::RegexMatch)
    print(io, "RegexMatch(")
    for k in 1:div(length(m.captured), 2)
        if k > 1
            print(io, ", ", k - 1, '=')
        end
        print(io, '"', m.seq[m.captured[2k-1]:m.captured[2k]-1], '"')
    end
    print(io, ')')
end

"""
    matched(match)

Return a matched pattern.
"""
function matched{S}(m::RegexMatch{S})
    return m.seq[m.captured[1]:m.captured[2]-1]
end

function matched{S}(m::Nullable{RegexMatch{S}})
    return matched(get(m))
end

"""
    captured(match)

Retrun a vector of captured patterns.
"""
function captured{S}(m::RegexMatch{S})
    return [m.captured[2k-1] != 0 && m.captured[2k] != 0 ?
            Nullable{S}(m.seq[m.captured[2k-1]:m.captured[2k]-1]) :
            Nullable{S}()
            for k in 2:div(length(m.captured), 2)]
end

function captured{S}(m::Nullable{RegexMatch{S}})
    return captured(get(m))
end

function checkeltype{T}(re::Regex{T}, seq::BioSequence)
    if eltype(seq) != T
        throw(ArgumentError("element type of sequence doesn't match with regex"))
    end
end

function Base.match{T}(re::Regex{T}, seq::BioSequence, start::Integer=1)
    checkeltype(re, seq)

    # use the first unambiguous symbol in the regular expression to find the
    # next starting position of pattern matching; this improves the performance
    if length(re.code) ≥ 2 && tag(re.code[2]) == BitsTag && count_ones(operand(re.code[2])) == 1
        firstsym = bits2sym(T, operand(re.code[2]))
    else
        firstsym = gap(T)
    end

    # a thread is `(<program counter>, <sequence's iterator state>)`
    threads = Stack{Tuple{Int,Int}}()
    captured = Vector{Int}(re.nsaves)
    s = start
    while true
        if firstsym != gap(T)
            s = findnext(seq, firstsym, s)
            if s == 0
                break
            end
        end
        empty!(threads)
        push!(threads, (1, s))
        fill!(captured, 0)
        if runmatch!(threads, captured, re, seq)
            return Nullable(RegexMatch(seq, captured))
        end
        _, s = next(seq, s)
        if done(seq, s)
            break
        end
    end
    return Nullable{RegexMatch{typeof(seq)}}()
end

function Base.search{T}(seq::BioSequence, re::Regex{T}, start::Integer=1)
    checkeltype(re, seq)
    m = Base.match(re, seq, start)
    if isnull(m)
        return 0:-1
    else
        x = get(m)
        return x.captured[1]:x.captured[2]-1
    end
end

immutable RegexMatchIterator{T,S}
    re::Regex{T}
    seq::S
    overlap::Bool

    function RegexMatchIterator(re::Regex{T}, seq::S, overlap::Bool)
        checkeltype(re, seq)
        return new(re, seq, overlap)
    end
end

function Base.iteratorsize{T,S}(::Type{RegexMatchIterator{T,S}})
    return Base.SizeUnknown()
end

function Base.eltype{T,S}(::Type{RegexMatchIterator{T,S}})
    return RegexMatch{S}
end

function Base.start(iter::RegexMatchIterator)
    threads = Stack{Tuple{Int,Int}}()
    captured = Vector{Int}(iter.re.nsaves)
    s = start(iter.seq)
    push!(threads, (1, s))
    fill!(captured, 0)
    return advance!(threads, captured, s, iter.re, iter.seq, iter.overlap)
end

function Base.done(iter::RegexMatchIterator, st)
    return isnull(st[1])
end

function Base.next(iter::RegexMatchIterator, st)
    item, threads, captured, s = st
    return get(item), advance!(threads, captured, s, iter.re, iter.seq, iter.overlap)
end

function advance!(threads, captured, s, re, seq, overlap)
    while true
        if runmatch!(threads, captured, re, seq)
            if !overlap
                empty!(threads)
                s = captured[2]
                if !done(seq, s)
                    push!(threads, (1, s))
                end
            end
            return Nullable(RegexMatch(seq, copy(captured))), threads, captured, s
        end
        if !isempty(seq)
            _, s = next(seq, s)
        end
        if done(seq, s)
            break
        end
        push!(threads, (1, s))
        fill!(captured, 0)
    end
    return Nullable{typeof(seq)}(), threads, captured, s
end

function Base.eachmatch{T}(re::Regex{T}, seq::BioSequence, overlap::Bool=true)
    checkeltype(re, seq)
    return RegexMatchIterator{T,typeof(seq)}(re, seq, overlap)
end

function Base.matchall{T}(re::Regex{T}, seq::BioSequence, overlap::Bool=true)
    # this will work on v0.5
    #   return map(matched, eachmatch(re, seq))
    ret = Vector{typeof(seq)}()
    for m in eachmatch(re, seq, overlap)
        push!(ret, matched(m))
    end
    return ret
end

function Base.ismatch{T}(re::Regex{T}, seq::BioSequence)
    return !isnull(Base.match(re, seq))
end

# simple stack
type Stack{T}
    top::Int
    data::Vector{T}

    Stack(sz::Int=0) = new(0, Vector{T}(sz))
end

function Base.isempty(stack::Stack)
    return stack.top == 0
end

@inline function Base.push!{T}(stack::Stack{T}, x::T)
    if stack.top + 1 > length(stack.data)
        push!(stack.data, x)
    else
        stack.data[stack.top+1] = x
    end
    stack.top += 1
    return stack
end

@inline function Base.pop!(stack::Stack)
    item = stack.data[stack.top]
    stack.top -= 1
    return item
end

@inline function Base.empty!(stack::Stack)
    stack.top = 0
    return stack
end

# run pattern maching using the current `threads` stack. Captured positions are
# stored in the preallocated `captured`. If match is found, this returns `true`
# immediately; otherwise returns `false`.
function runmatch!(threads::Stack{Tuple{Int,Int}},
                   captured::Vector{Int},
                   re::Regex, seq::BioSequence)
    while !isempty(threads)
        pc::Int, s = pop!(threads)
        while true
            op = re.code[pc]
            t = tag(op)
            if t == BitsTag
                if done(seq, s)
                    break
                end
                sym, s = next(seq, s)
                if Seq.compatbits(sym) & operand(op) != 0
                    pc += 1
                else
                    break
                end
            elseif t == JumpTag
                pc = operand(op)
            elseif t == PushTag
                push!(threads, (convert(Int, operand(op)), s))
                pc += 1
            elseif t == ForkTag
                push!(threads, (pc + 1, s))
                pc = operand(op)
            elseif t == SaveTag
                captured[operand(op)] = s
                pc += 1
            elseif t == HeadTag
                if s == start(seq)
                    pc += 1
                else
                    break
                end
            elseif t == LastTag
                if s == endof(seq) + 1
                    pc += 1
                else
                    break
                end
            elseif t == MatchTag
                return true
            end
        end
    end
    return false
end

end  # module RE

# exported from Bio.Seq
const matched = RE.matched
const captured = RE.captured
