# Precompile
# ==========

# Parser
# ------

precompile(Base.open, (String, Type{Seq.FASTA},))
precompile(Base.open, (String, Type{Seq.FASTQ},))
precompile(Base.open, (String, Type{Intervals.BED},))
precompile(Base.read, (Seq.FASTAParser{Seq.BioSequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.DNASequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.RNASequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.AminoAcidSequence},))
precompile(Base.read, (Seq.FASTQParser{Seq.DNASequence},))
precompile(Base.read, (Intervals.BEDParser,))
