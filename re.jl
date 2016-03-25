# Regular Expression
# ==================

module RE

using Bio.Seq


# Syntax tree
# -----------

type SyntaxTree
    head::Symbol
    args::Vector{Any}

    function SyntaxTree(head, args)
        @assert head ∈ (
            :*, :+, :?, :|,
            :range, :set, :compset, :sym, :bits,
            :capture, :concat, :head, :last)
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

const symbols = Dict(
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

parse{T}(::Type{T}, pat::AbstractString) = parserec(T, pat, start(pat))[1]

function parserec{T}(::Type{T}, pat, s)
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
            push!(args, expr(:?, [arg]))
        elseif c == '{'
            @check !isempty(args) ArgumentError("unexpected '{'")
            rng, s = parserange(pat , s)
            arg = pop!(args)
            push!(args, expr(:range, [rng, arg]))
        elseif c == '|'
            arg1 = expr(:concat, args)
            arg2, s = parserec(T, pat, s)
            args = []
            push!(args, expr(:|, [arg1, arg2]))
        elseif c == '['
            setexpr, s = parseset(T, pat, s)
            push!(args, setexpr)
        elseif c == '('
            arg, s = parserec(T, pat, s)
            push!(args, expr(:capture, [arg]))
        elseif c == ')'
            return expr(:concat, args), s
        elseif c == '^'
            push!(args, expr(:head, []))
        elseif c == '$'
            push!(args, expr(:last, []))
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

function parseset{T}(::Type{T}, pat, s)
    set = T[]
    complement = false
    while !done(pat, s)
        c, s = next(pat, s)
        if c ∈ symbols[T]
            push!(set, convert(T, c))
        elseif c == '^'
            complement = true
        elseif c == ']'
            break
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    return expr(complement ? :compset : :set, set), s
end

# DNA/RNA nucleotides
const sym2bits_nuc = UInt32[
    0b0001,
    0b0010,
    0b0100,
    0b1000,
    0b0011,
    0b0101,
    0b1001,
    0b0110,
    0b1010,
    0b1100,
    0b0111,
    0b1011,
    0b1101,
    0b1110,
    0b1111]
sym2bits(nt::DNANucleotide) = sym2bits_nuc[reinterpret(Int8, nt)+1]
sym2bits(nt::RNANucleotide) = sym2bits_nuc[reinterpret(Int8, nt)+1]
mask(::Type{DNANucleotide}) = (UInt32(1) << 4) - one(UInt32)
mask(::Type{RNANucleotide}) = (UInt32(1) << 4) - one(UInt32)

const sym2bits_aa = UInt32[
    0b0000000000000000000001,  # A
    0b0000000000000000000010,  # R
    0b0000000000000000000100,  # N
    0b0000000000000000001000,  # D
    0b0000000000000000010000,  # C
    0b0000000000000000100000,  # Q
    0b0000000000000001000000,  # E
    0b0000000000000010000000,  # G
    0b0000000000000100000000,  # H
    0b0000000000001000000000,  # I
    0b0000000000010000000000,  # L
    0b0000000000100000000000,  # K
    0b0000000001000000000000,  # M
    0b0000000010000000000000,  # F
    0b0000000100000000000000,  # P
    0b0000001000000000000000,  # S
    0b0000010000000000000000,  # T
    0b0000100000000000000000,  # W
    0b0001000000000000000000,  # Y
    0b0010000000000000000000,  # V
    0b0100000000000000000000,  # O
    0b1000000000000000000000,  # U
    0b0000000000000000001100,  # B
    0b0000000000011000000000,  # J
    0b0000000000000001100000,  # Z
    0b1111111111111111111111,  # X
]
sym2bits(aa::AminoAcid) = sym2bits_aa[Int8(aa)+1]
mask(::Type{AminoAcid}) = (UInt32(1) << 22) - one(UInt32)

function desugar{T}(::Type{T}, tree::SyntaxTree)
    head = tree.head
    args = tree.args
    if head == :+
        # e+ => ee*
        head = :concat
        args = [args[1], expr(:*, [args[1]])]
    elseif head == :?
        # e? => e|
        head = :|
        args = [args[1], expr(:concat, [])]
    elseif head == :sym
        head = :bits
        args = [sym2bits(args[1])]
    elseif head == :set
        bits = UInt32(0)
        for arg in args
            bits |= sym2bits(arg)
        end
        head = :bits
        args = [bits]
    elseif head == :compset
        bits = UInt32(0)
        for arg in args
            bits |= sym2bits(arg)
        end
        head = :bits
        args = [~bits & mask(T)]
    elseif head == :range
        rng = args[1]
        pat = args[2]
        if isa(rng, Int)
            # e{m} => eee...e
            #         |<-m->|
            head = :concat
            args = replicate(pat, rng)
        elseif isa(rng, Tuple{Int})
            # e{m,} => eee...ee*
            #          |<-m->|
            head = :concat
            args = replicate(pat, rng[1])
            push!(args, expr(:*, [pat]))
        elseif isa(rng, Tuple{Int,Int})
            # e{m,n} => eee...e|eee...ee|...|eee...eee
            #           |<-m->|              |<- n ->|
            m, n = rng
            if m == n
                head = :concat
                args = replicate(pat, m)
            else
                @assert m < n
                tree = replicate(pat, n)
                for k in n-1:-1:m
                    tree = expr(:|, [replicate(pat, k), tree])
                end
                head = :|
                args = tree.args
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

desugar{T}(::Type{T}, atom) = atom

replicate(x, n) = collect(repeated(x, n))


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
tag(op::Op) = reinterpret(UInt32, op) & ~operand_mask
operand(op::Op) = reinterpret(UInt32, op) & operand_mask

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
        push!(code, bits(args[1]))
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

immutable Regex{T}
    pat::ASCIIString  # regular expression pattern (for printing)
    code::Vector{Op}  # compiled code
    nsaves::Int       # the number of `save` operations in `code`

    function Regex(pat::AbstractString)
        code = compile(desugar(T, parse(T, pat)))
        nsaves = 0
        for op in code
            nsaves += tag(op) == SaveTag
        end
        @assert iseven(nsaves)
        return new(pat, code, nsaves)
    end
end

function Base.show(io::IO, re::Regex)
    print(io, summary(re), "(\"", re.pat, "\")")
end

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

function matched(m::RegexMatch)
    return m.seq[m.captured[1]:m.captured[2]-1]
end

function captured(m::RegexMatch)
    return [m.seq[m.captured[2k-1]:m.captured[2k]-1]
            for k in 2:div(length(m.captured), 2)]
end

# simple stack
type Stack{T}
    top::Int
    data::Vector{T}

    Stack(sz::Int=0) = new(0, Vector{T}(sz))
end

Base.isempty(stack::Stack) = stack.top == 0

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

function match{T}(re::Regex{T}, seq::BioSequence)
    if eltype(seq) != T
        throw(ArgumentError("element type of sequence doesn't match with regex"))
    end
    ss = start(seq)
    # a thread is `(<program counter>, <sequence's iterator state>)`
    threads = Stack{Tuple{Int,Int}}()
    captured = zeros(typeof(ss), re.nsaves)
    while !done(seq, ss)
        push!(threads, (1, ss))
        while !isempty(threads)
            pc, s = pop!(threads)
            while true
                op = re.code[pc]
                t = tag(op)
                if t == BitsTag
                    if done(seq, s)
                        break
                    end
                    sym, s = next(seq, s)
                    if sym2bits(sym) & operand(op) != 0
                        pc += 1
                    else
                        break
                    end
                elseif t == JumpTag
                    pc = convert(Int, operand(op))
                elseif t == PushTag
                    push!(threads, (convert(Int, operand(op)), s))
                    pc += 1
                elseif t == ForkTag
                    push!(threads, (pc + 1, s))
                    pc = convert(Int, operand(op))
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
                    return Nullable(RegexMatch(seq, captured))
                end
            end
        end
        empty!(threads)
        _, ss = next(seq, ss)
    end
    return Nullable{RegexMatch{eltype(seq)}}()
end

ismatch{T}(re::Regex{T}, seq::BioSequence) = !isnull(match(re, seq))

# inline quick tests
using Base.Test
@test  ismatch(Regex{DNANucleotide}("A"), dna"AA")
@test  ismatch(Regex{DNANucleotide}("A+"), dna"AA")
@test !ismatch(Regex{DNANucleotide}("A+"), dna"CC")
@test  ismatch(Regex{DNANucleotide}("A+C+"), dna"AAC")
@test !ismatch(Regex{DNANucleotide}("A+C+"), dna"AA")
@test  ismatch(Regex{DNANucleotide}("T(A[AG]|GA)"), dna"TGA")

@test  matched(get(match(Regex{DNANucleotide}("A(C+)"), dna"ACCC"))) == dna"ACCC"
@test captured(get(match(Regex{DNANucleotide}("A(C+)"), dna"ACCC"))) == [dna"CCC"]
@test captured(get(match(Regex{DNANucleotide}("(A)(C+)"), dna"ACCC"))) == [dna"A", dna"CCC"]

end  # module RE
