# Genetic Code
# ============
#
# Genetic code table and translator from RNA to amino acid sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# A genetic code is a table mapping RNA 3-mers (i.e. RNAKmer{3}) to AminoAcids.
"Type representing a Genetic Code"
immutable GeneticCode <: Associative{RNACodon, AminoAcid}
    name::Compat.String
    tbl::Vector{AminoAcid}
end


# Basic Functions
# ---------------

function Base.getindex(code::GeneticCode, idx::Union{DNACodon,RNACodon})
    return code.tbl[convert(UInt64, idx) + 1]
end

function Base.setindex!(code::GeneticCode, aa::AminoAcid, idx::Union{DNACodon,RNACodon})
    return setindex!(code.tbl, aa, convert(UInt64, idx) + 1)
end

Base.copy(code::GeneticCode) = GeneticCode(copy(code.name), copy(code.tbl))
Base.length(code::GeneticCode) = 64

Base.showcompact(io::IO, code::GeneticCode) = print(io, code.name)

function Base.show(io::IO, code::GeneticCode)
    print(io, code.name)
    rna = rna"ACGU"
    for x in rna, y in rna
        println(io)
        print(io, "  ")
        for z in rna
            codon = Kmer(x, y, z)
            aa = code[codon]
            print(io, codon, ": ", aa)
            if z != RNA_U
                print(io, "    ")
            end
        end
    end
end


# Iterating through genetic code
# ------------------------------

Base.start(code::GeneticCode) = UInt64(0)

function Base.next(code::GeneticCode, x::UInt64)
    c = convert(RNACodon, x)
    return (c, code[c]), x + 1
end

Base.done(code::GeneticCode, x::UInt64) = x > UInt64(0b111111)


# Default genetic codes
# ---------------------

immutable TransTables
    tables::Dict{Int,GeneticCode}
    bindings::Dict{Int,Symbol}
    function TransTables()
        return new(Dict(), Dict())
    end
end

Base.getindex(trans::TransTables, key::Integer) = trans.tables[Int(key)]

function Base.show(io::IO, trans::TransTables)
    print(io, "Translation Tables:")
    ids = sort(collect(keys(trans.tables)))
    for id in ids
        println(io)
        print(io, lpad(id, 3), ". ")
        showcompact(io, trans.tables[id])
        if haskey(trans.bindings, id)
            print(io, " (", trans.bindings[id], ")")
        end
    end
end

"""
Genetic code list of NCBI.

The standard genetic code is `ncbi_trans_table[1]` and others can be shown by
`show(ncbi_trans_table)`.
For more details, consult the next link:
http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes.
"""
const ncbi_trans_table = TransTables()

macro register_ncbi_gencode(id, bind, tbl)
    quote
        gencode = parse_gencode($tbl)
        const $(esc(bind)) = gencode
        ncbi_trans_table.tables[$id] = gencode
        ncbi_trans_table.bindings[$id] = Symbol($(string(bind)))
    end
end

function parse_gencode(s)
    name, _, aas, _, base1, base2, base3 = split(chomp(s), '\n')
    name = split(name, ' ', limit=2)[2]  # drop number
    codes = GeneticCode(name, Vector{AminoAcid}(4^3))
    @assert length(aas) == 73
    for i in 10:73
        aa = AminoAcid(aas[i])
        b1 = DNANucleotide(base1[i])
        b2 = DNANucleotide(base2[i])
        b3 = DNANucleotide(base3[i])
        codon = DNACodon(Kmer(b1, b2, b3))
        codes[codon] = aa
    end
    return codes
end

# Genetic codes translation tables are taken from the NCBI taxonomy database.

@register_ncbi_gencode 1 standard_genetic_code """
1. The Standard Code

  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 2 vertebrate_mitochondrial_genetic_code """
2. The Vertebrate Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
Starts = --------------------------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 3 yeast_mitochondrial_genetic_code """
3. The Yeast Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ----------------------------------MM----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 4 mold_mitochondrial_genetic_code """
4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = --MM---------------M------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 5 invertebrate_mitochondrial_genetic_code """
5. The Invertebrate Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
Starts = ---M----------------------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 6 ciliate_nuclear_genetic_code """
6. The Ciliate, Dasycladacean and Hexamita Nuclear Code

  AAs  = FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 9 echinoderm_mitochondrial_genetic_code """
9. The Echinoderm and Flatworm Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 10 euplotid_nuclear_genetic_code """
10. The Euplotid Nuclear Code

  AAs  = FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 11 bacterial_plastid_genetic_code """
11. The Bacterial, Archaeal and Plant Plastid Code

  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 12 alternative_yeast_nuclear_genetic_code """
12. The Alternative Yeast Nuclear Code

  AAs  = FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -------------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 13 ascidian_mitochondrial_genetic_code """
13. The Ascidian Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
Starts = ---M------------------------------MM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 14 alternative_flatworm_mitochondrial_genetic_code """
14. The Alternative Flatworm Mitochondrial Code

  AAs  = FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 16 chlorophycean_mitochondrial_genetic_code """
16. Chlorophycean Mitochondrial Code

  AAs  = FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 21 trematode_mitochondrial_genetic_code """
21. Trematode Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 22 scenedesmus_obliquus_mitochondrial_genetic_code """
22. Scenedesmus obliquus Mitochondrial Code

  AAs  = FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 23 thraustochytrium_mitochondrial_genetic_code """
23. Thraustochytrium Mitochondrial Code

  AAs  = FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = --------------------------------M--M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 24 pterobrachia_mitochondrial_genetic_code """
24. Pterobranchia Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 25 candidate_division_sr1_genetic_code """
25. Candidate Division SR1 and Gracilibacteria Code

  AAs  = FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M-------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""


# Translation
# -----------

"""
    translate(rna_seq, code=standard_genetic_code, allow_ambiguous_codons=true)

Translate an `RNASequence` to an `AminoAcidSequence`.

Translation uses genetic code `code` to map codons to amino acids. See
`ncbi_trans_table` for available genetic codes.
If codons in the given RNA sequence cannot determine a unique amino acid, they
will be translated to `AA_X` if `allow_ambiguous_codons` is `true` and otherwise
result in an error.
"""
function translate(seq::RNASequence;
                   code::GeneticCode=standard_genetic_code,
                   allow_ambiguous_codons::Bool=true)
    return translate(seq, code, allow_ambiguous_codons)
end

function translate(seq::RNASequence, code::GeneticCode, allow_ambiguous_codons::Bool)
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
            aaseq[j] = code[Kmer(x, y, z)]
        end
        i += 3
        j += 1
    end
    return aaseq
end

function try_translate_ambiguous_codon(code::GeneticCode,
                                       x::RNANucleotide,
                                       y::RNANucleotide,
                                       z::RNANucleotide)
    if !isambiguous(x) && !isambiguous(y)
        # try to translate a codon `(x, y, RNA_N)`
        aa_a = code[Kmer(x, y, RNA_A)]
        aa_c = code[Kmer(x, y, RNA_C)]
        aa_g = code[Kmer(x, y, RNA_G)]
        aa_u = code[Kmer(x, y, RNA_U)]
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
