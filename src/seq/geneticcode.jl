# Genetic Code
# ============

# A genetic code is a table mapping RNA 3-mers (i.e. RNAKmer{3}) to AminoAcids.
"Type representing a Genetic Code"
immutable GeneticCode <: Associative{RNAKmer{3}, AminoAcid}
    tbl::Vector{AminoAcid}
end


# Basic Functions
# ---------------

@inline getindex(code::GeneticCode, idx::RNAKmer{3}) = code.tbl[convert(Uint64, idx) + 1]
setindex!(code::GeneticCode, aa::AminoAcid, idx::RNAKmer{3}) = (code.tbl[convert(Uint64, idx) + 1] = aa)
copy(code::GeneticCode) = GeneticCode(copy(code.tbl))
length(code::GeneticCode) = 64


# Iterating through genetic code
# ------------------------------

start(code::GeneticCode) = @compat uint64(0)

function next(code::GeneticCode, x::Uint64)
    c = convert(Codon, x)
    return ((c, code[c]), (x + 1))
end

done(code::GeneticCode, x::Uint64) = (x > (@compat uint64(0b111111)))


# Default genetic codes
# ---------------------
# All of these taken from:
# http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1

# Standard Code
const standard_genetic_code = GeneticCode([
    AA_K, AA_Q, AA_E, AA_INVALID, AA_T, AA_P, AA_A, AA_S, AA_R, AA_R, AA_G,
    AA_INVALID, AA_I, AA_L, AA_V, AA_L, AA_N, AA_H, AA_D, AA_Y, AA_T, AA_P,
    AA_A, AA_S, AA_S, AA_R, AA_G, AA_C, AA_I, AA_L, AA_V, AA_F, AA_K, AA_Q,
    AA_E, AA_INVALID, AA_T, AA_P, AA_A, AA_S, AA_R, AA_R, AA_G, AA_W, AA_M,
    AA_L, AA_V, AA_L, AA_N, AA_H, AA_D, AA_Y, AA_T, AA_P, AA_A, AA_S, AA_S,
    AA_R, AA_G, AA_C, AA_I, AA_L, AA_V, AA_F ])

# Vertebrate mithocondrial genetic code
const vertebrate_mitochondrial_genetic_code = copy(standard_genetic_code)
vertebrate_mitochondrial_genetic_code[rnakmer("AGA")] = AA_INVALID
vertebrate_mitochondrial_genetic_code[rnakmer("AGG")] = AA_INVALID
vertebrate_mitochondrial_genetic_code[rnakmer("AUA")] = AA_M
vertebrate_mitochondrial_genetic_code[rnakmer("UGA")] = AA_W

# Yeast mithocondrial genetic code
const yeast_mitochondrial_genetic_code = copy(standard_genetic_code)
yeast_mitochondrial_genetic_code[rnakmer("AUA")] = AA_M
yeast_mitochondrial_genetic_code[rnakmer("CUU")] = AA_T
yeast_mitochondrial_genetic_code[rnakmer("CUC")] = AA_T
yeast_mitochondrial_genetic_code[rnakmer("CUA")] = AA_T
yeast_mitochondrial_genetic_code[rnakmer("CUG")] = AA_T
yeast_mitochondrial_genetic_code[rnakmer("UGA")] = AA_W
yeast_mitochondrial_genetic_code[rnakmer("CGA")] = AA_INVALID
yeast_mitochondrial_genetic_code[rnakmer("CGC")] = AA_INVALID

# Mold mithocondrial genetic code
const mold_mitochondrial_genetic_code = copy(standard_genetic_code)
mold_mitochondrial_genetic_code[rnakmer("UGA")] = AA_W

# Invertebrate mithocondrial genetic code
const invertebrate_mitochondrial_genetic_code = copy(standard_genetic_code)
invertebrate_mitochondrial_genetic_code[rnakmer("AGA")] = AA_S
invertebrate_mitochondrial_genetic_code[rnakmer("AGG")] = AA_S
invertebrate_mitochondrial_genetic_code[rnakmer("AUA")] = AA_M
invertebrate_mitochondrial_genetic_code[rnakmer("AGA")] = AA_W

# Ciliate nuclear genetic code
const ciliate_nuclear_genetic_code = copy(standard_genetic_code)
ciliate_nuclear_genetic_code[rnakmer("UAA")] = AA_Q
ciliate_nuclear_genetic_code[rnakmer("UAG")] = AA_Q

# Echinoderm mithocondrial genetic code
const echinoderm_mitochondrial_genetic_code = copy(standard_genetic_code)
echinoderm_mitochondrial_genetic_code[rnakmer("AAA")] = AA_N
echinoderm_mitochondrial_genetic_code[rnakmer("AGA")] = AA_S
echinoderm_mitochondrial_genetic_code[rnakmer("AGG")] = AA_S
echinoderm_mitochondrial_genetic_code[rnakmer("UGA")] = AA_W

# Euplotid nuclear genetic code
const euplotid_nuclear_genetic_code = copy(standard_genetic_code)
euplotid_nuclear_genetic_code[rnakmer("UGA")] = AA_C

# Bacterial plastid genetic code
const bacterial_plastid_genetic_code = copy(standard_genetic_code)

# Alternative yeast nuclear genetic code
const alternative_yeast_nuclear_genetic_code = copy(standard_genetic_code)
alternative_yeast_nuclear_genetic_code[rnakmer("CUG")] = AA_S

# Ascidian mithocondrial genetic code
const ascidian_mitochondrial_genetic_code = copy(standard_genetic_code)
ascidian_mitochondrial_genetic_code[rnakmer("AGA")] = AA_G
ascidian_mitochondrial_genetic_code[rnakmer("AGG")] = AA_G
ascidian_mitochondrial_genetic_code[rnakmer("AUA")] = AA_M
ascidian_mitochondrial_genetic_code[rnakmer("UGA")] = AA_W

# Alternative flatworm mithocondrial genetic code
const alternative_flatworm_mitochondrial_genetic_code = copy(standard_genetic_code)
alternative_flatworm_mitochondrial_genetic_code[rnakmer("AAA")] = AA_N
alternative_flatworm_mitochondrial_genetic_code[rnakmer("AGA")] = AA_S
alternative_flatworm_mitochondrial_genetic_code[rnakmer("AGG")] = AA_S
alternative_flatworm_mitochondrial_genetic_code[rnakmer("UAA")] = AA_Y
alternative_flatworm_mitochondrial_genetic_code[rnakmer("UGA")] = AA_W

# Clorophycean mithocondrial genetic code
const chlorophycean_mitochondrial_genetic_code = copy(standard_genetic_code)
chlorophycean_mitochondrial_genetic_code[rnakmer("UAG")] = AA_L

# Trematode mithocondrial genetic code
const trematode_mitochondrial_genetic_code = copy(standard_genetic_code)
trematode_mitochondrial_genetic_code[rnakmer("UGA")] = AA_W
trematode_mitochondrial_genetic_code[rnakmer("AUA")] = AA_M
trematode_mitochondrial_genetic_code[rnakmer("AGA")] = AA_S
trematode_mitochondrial_genetic_code[rnakmer("AGG")] = AA_S
trematode_mitochondrial_genetic_code[rnakmer("AAA")] = AA_N

# Scenedesmus obliquus mithocondrial genetic code
const scenedesmus_obliquus_mitochondrial_genetic_code = copy(standard_genetic_code)
scenedesmus_obliquus_mitochondrial_genetic_code[rnakmer("UCA")] = AA_INVALID
scenedesmus_obliquus_mitochondrial_genetic_code[rnakmer("UAG")] = AA_L

# Vertebrate mithocondrial genetic code
const thraustochytrium_mitochondrial_genetic_code = copy(standard_genetic_code)
thraustochytrium_mitochondrial_genetic_code[rnakmer("UUA")] = AA_INVALID

# Pterobrachia mithocondrial genetic code
const pterobrachia_mitochondrial_genetic_code = copy(standard_genetic_code)
pterobrachia_mitochondrial_genetic_code[rnakmer("AGA")] = AA_S
pterobrachia_mitochondrial_genetic_code[rnakmer("AGG")] = AA_K
pterobrachia_mitochondrial_genetic_code[rnakmer("UGA")] = AA_W

# Candidate division sr1 genetic code
const candidate_division_sr1_genetic_code = copy(standard_genetic_code)
candidate_division_sr1_genetic_code[rnakmer("UGA")] = AA_INVALID

# Genetic codes indexed as in NCBI's trans_table
const ncbi_trans_table = GeneticCode[
    standard_genetic_code,
    vertebrate_mitochondrial_genetic_code,
    yeast_mitochondrial_genetic_code,
    mold_mitochondrial_genetic_code,
    invertebrate_mitochondrial_genetic_code,
    ciliate_nuclear_genetic_code,
    echinoderm_mitochondrial_genetic_code,
    euplotid_nuclear_genetic_code,
    bacterial_plastid_genetic_code,
    alternative_yeast_nuclear_genetic_code,
    ascidian_mitochondrial_genetic_code,
    alternative_flatworm_mitochondrial_genetic_code,
    chlorophycean_mitochondrial_genetic_code,
    trematode_mitochondrial_genetic_code,
    scenedesmus_obliquus_mitochondrial_genetic_code,
    thraustochytrium_mitochondrial_genetic_code,
    pterobrachia_mitochondrial_genetic_code,
    candidate_division_sr1_genetic_code ]


# Translation
# -----------

# Try to translate a codon (u, v, RNA_N). Return the corresponding amino acid,
# if it can be unambiguously translated, and AA_INVALID if it cannot be.
function translate_ambiguous_codon(code::GeneticCode, u::RNANucleotide,
                                   v::RNANucleotide)
    aa_a = code[kmer(u, v, RNA_A)]
    aa_c = code[kmer(u, v, RNA_C)]
    aa_g = code[kmer(u, v, RNA_G)]
    aa_u = code[kmer(u, v, RNA_U)]

    if aa_a == aa_c == aa_g == aa_u
        return aa_a
    else
        error("Codon $(u)$(v)N cannot be unambiguously translated.")
    end
end


"Convert an RNASequence to an AminoAcidSequence"
function translate(seq::RNASequence, code::GeneticCode=standard_genetic_code,
                   allow_ambiguous_codons::Bool=false)
    aaseqlen, r = divrem(length(seq), 3)
    if r != 0
        error("RNASequence length is not divisible by three. Cannot translate.")
    end

    # try to translate codons whose third nucleotide is N
    aaseq = Array(AminoAcid, aaseqlen)

    for i in npositions(seq)
        d, r = divrem(i - 1, 3)

        if r != 2
            if allow_ambiguous_codons
                aaseq[d + 1] = AA_X
            else
                if r == 0
                    codon_str = "N$(seq[i+1])$(seq[i+2])"
                else
                    codon_str = "$(seq[i-1])N$(seq[i+1])"
                end
                error("Codon $(codon_str) cannot be unambiguously translated.")
            end
        else
            aa = translate_ambiguous_codon(code, seq[i-2], seq[i-1])
            aaseq[d + 1] = aa
        end
    end
    @inbounds for (i, codon) in each(RNAKmer{3}, seq, 3)
        aa = code[codon]
        if aa == AA_INVALID
            error("Cannot translate stop codons.")
        end
        j = div(i - 1, 3) + 1
        aaseq[j] = code[codon]
    end

    return AminoAcidSequence(aaseq, 1:aaseqlen)
end
