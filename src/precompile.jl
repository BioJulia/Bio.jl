# Precompile
# ==========

# Parser
# ------

if VERSION < v"0.5-"
    precompile(Base.open, (ASCIIString, Type{Seq.FASTA},))
    precompile(Base.open, (ASCIIString, Type{Seq.FASTQ}, Type{Seq.QualityEncoding},))
    precompile(Base.open, (ASCIIString, Type{Intervals.BED},))
else
    precompile(Base.open, (String, Type{Seq.FASTA},))
    precompile(Base.open, (String, Type{Seq.FASTQ}, Type{Seq.QualityEncoding},))
    precompile(Base.open, (String, Type{Intervals.BED},))
end
precompile(Base.read, (Seq.FASTAParser{Seq.BioSequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.DNASequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.RNASequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.AminoAcidSequence},))
precompile(Base.read, (Seq.FASTQParser{Seq.DNASequence},))
precompile(Base.read, (Intervals.BEDParser,))
