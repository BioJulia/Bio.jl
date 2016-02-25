# Multinominal naive Bayes classifier to infer alphabet from sequence
# See https://en.wikipedia.org/wiki/Naive_Bayes_classifier for details

# A/C/G/T/U/N/Other
const n_features = 7
const alphabets = [DNAAlphabet, RNAAlphabet, AminoAcidAlphabet]
const alphabet_type = Dict(
    DNAAlphabet => DNASequence,
    RNAAlphabet => RNASequence,
    AminoAcidAlphabet => AminoAcidSequence)

# logarithmic likelihood and prior for each character given an alphabet
# NOTE: Likelihood was trained with the train function.
const loglikelihood = [
    2.856e-01 2.066e-01 7.765e-02;
    1.979e-01 2.799e-01 7.355e-03;
    2.141e-01 3.066e-01 9.444e-02;
    3.024e-01 6.664e-05 6.611e-02;
    1.999e-06 2.066e-01 1.049e-05;
    1.999e-06 6.664e-05 2.834e-02;
    1.999e-06 6.664e-05 7.261e-01;
]
#const logprior = log(ones(length(alphabets)) / length(alphabets))
const logprior = log([0.45, 0.1, 0.45])

function predict(seq::Vector{UInt8}, start, stop)
    # count characters
    a = c = g = t = u = n = other = 0
    #for x in seq
    for i in start:stop
        x = seq[i]
        if x == 0x41 || x == 0x61
            a += 1
        elseif x == 0x43 || x == 0x63
            c += 1
        elseif x == 0x47 || x == 0x67
            g += 1
        elseif x == 0x54 || x == 0x74
            t += 1
        elseif x == 0x55 || x == 0x75
            u += 1
        elseif x == 0x4e || x == 0x6e
            n += 1
        else
            other += 1
        end
    end
    # compute logarithmic posterior for each alphabet
    K = length(logprior)
    logposterior = copy(logprior)
    for k in 1:K
        logposterior[k] +=
            loglikelihood[1,k] * a +
            loglikelihood[2,k] * c +
            loglikelihood[3,k] * g +
            loglikelihood[4,k] * t +
            loglikelihood[5,k] * u +
            loglikelihood[6,k] * n +
            loglikelihood[7,k] * other
    end
    return alphabets[indmax(logposterior)]
end

predict(seq::ASCIIString) = predict(seq.data)

function train(seqs::Vector{ASCIIString}, labels::Vector{Int};
               pseudo_count::Float64=0.01)
    @assert length(seqs) == length(labels)
    K = length(alphabets)
    counts = [ones(n_features) * pseudo_count for _ in 1:K]
    for n in 1:endof(seqs)
        for x in seqs[n]
            x = lowercase(x)
            f = x == 'a' ? 1 :
                x == 'c' ? 2 :
                x == 'g' ? 3 :
                x == 't' ? 4 :
                x == 'u' ? 5 :
                x == 'n' ? 6 : 7
            counts[labels[n]][f] += 1
        end
    end
    for k in 1:K
        scale!(counts[k], 1 / sum(counts[k]))
        @assert sum(counts[k]) ≈ 1
    end
    return hcat(counts...)
end

function dotrain(filepath)
    open(filepath) do f
        seqs = ASCIIString[]
        labels = Int[]
        seq = ""
        label = 0
        for line in eachline(f)
            if line[1] == '>'
                if !isempty(seq)
                    push!(seqs, seq)
                    push!(labels, label)
                    seq = ""
                end
                desc = chomp(line[2:end])
                label = desc == "dna"       ? 1 :
                        desc == "rna"       ? 2 :
                        desc == "aminoacid" ? 3 :
                        error("unknown label: $desc")
            else
                seq *= chomp(line)
            end
        end
        if !isempty(seq)
            push!(seqs, seq)
            push!(labels, label)
        end
        P = train(seqs, labels)
        println("const loglikelihood = [")
        for i in 1:size(P, 1)
            print("    ")
            for j in 1:size(P, 2)
                @printf "%.3e" P[i,j]
                if j == size(P, 2)
                    print(";\n")
                else
                    print(" ")
                end
            end
        end
        println("]")
    end
end
