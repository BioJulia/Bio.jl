# Genetic Code
# ============

# A genetic code is a table mapping RNA 3-mers (i.e. RNAKmer{3}) to AminoAcids.
"Type representing a Genetic Code"
immutable GeneticCode <: Associative{RNAKmer{3}, AminoAcid}
    tbl::Vector{AminoAcid}
end


# Basic Functions
# ---------------

@inline Base.getindex(code::GeneticCode, idx::RNAKmer{3}) = code.tbl[convert(UInt64, idx) + 1]
Base.setindex!(code::GeneticCode, aa::AminoAcid, idx::RNAKmer{3}) = (code.tbl[convert(UInt64, idx) + 1] = aa)
Base.copy(code::GeneticCode) = GeneticCode(copy(code.tbl))
Base.length(code::GeneticCode) = 64


# Iterating through genetic code
# ------------------------------

Base.start(code::GeneticCode) = UInt64(0)

function Base.next(code::GeneticCode, x::UInt64)
    c = convert(Codon, x)
    return ((c, code[c]), (x + 1))
end

Base.done(code::GeneticCode, x::UInt64) = (x > UInt64(0b111111))


# Default genetic codes
# ---------------------
# All of these taken from:
# http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1

function parse_genetic_codes(s)
    aas, _, base1, base2, base3 = split(chomp(s), '\n')
    codes = GeneticCode(Vector{AminoAcid}(4^3))
    @assert length(aas) == 73
    for i in 10:73
        aa = AminoAcid(aas[i])
        b1 = DNANucleotide(base1[i])
        b2 = DNANucleotide(base2[i])
        b3 = DNANucleotide(base3[i])
        codon = Codon(kmer(b1, b2, b3))
        codes[codon] = aa
    end
    return codes
end

# Standard Code
const standard_genetic_code = parse_genetic_codes("""
  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
""")

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

function try_translate_ambiguous_codon(code::GeneticCode,
                                       x::RNANucleotide,
                                       y::RNANucleotide,
                                       z::RNANucleotide)
    if !isambiguous(x) && !isambiguous(y)
        # try to translate a codon `(x, y, RNA_N)`
        aa_a = code[kmer(x, y, RNA_A)]
        aa_c = code[kmer(x, y, RNA_C)]
        aa_g = code[kmer(x, y, RNA_G)]
        aa_u = code[kmer(x, y, RNA_U)]
        if aa_a == aa_c == aa_g == aa_u
            return Nullable{AminoAcid}(aa_a)
        end
    end

    found = Nullable{AminoAcid}()
    for (codon, aa) in code
        if (iscompatible(x, codon[1]) &&
            iscompatible(y, codon[2]) &&
            iscompatible(z, codon[3]))
            if isnull(found)
                found = Nullable(aa)
            elseif aa != get(found)
                return Nullable{AminoAcid}()
            end
        end
    end
    return found
end


"""
Translate an `RNASequence` to an `AminoAcidSequence`.

### Arguments
  * `seq`: RNA sequence to translate.
  * `code`: Genetic code to use (default is the standard genetic code).
  * `allow_ambiguous_codons`: True if ambiguous codons should be allowed and
      translated to `AA_X`. If false, they will throw an error. (default is true)

### Returns
A translated `AminoAcidSequence`
"""
function translate(seq::RNASequence,
                   code::GeneticCode=standard_genetic_code,
                   allow_ambiguous_codons::Bool=true)
    aaseqlen, r = divrem(length(seq), 3)
    if r != 0
        error("RNASequence length is not divisible by three. Cannot translate.")
    end

    aaseq = AminoAcidSequence(aaseqlen)
    i = j = 1
    while i ≤ endof(seq) - 2
        x = seq[i]
        y = seq[i+1]
        z = seq[i+2]
        if isambiguous(x) || isambiguous(y) || isambiguous(z)
            aa = try_translate_ambiguous_codon(code, x, y, z)
            if isnull(aa)
                if allow_ambiguous_codons
                    aaseq[j] = AA_X
                else
                    error("codon ", x, y, z, " cannot be unambiguously translated")
                end
            else
                aaseq[j] = get(aa)
            end
        else
            aaseq[j] = code[kmer(x, y, z)]
        end
        i += 3
        j += 1
    end
    return aaseq
end

function translate(seq::RNASequence;
                   code::GeneticCode=standard_genetic_code,
                   allow_ambiguous_codons::Bool=true)
    return translate(seq, code, allow_ambiguous_codons)
end
