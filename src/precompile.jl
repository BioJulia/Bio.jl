# Precompile
# ==========

# Parser
# ------

precompile(Base.open, (ASCIIString, Type{Seq.FASTA},))
precompile(Base.read, (Seq.FASTAParser{Seq.BioSequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.DNASequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.RNASequence},))
precompile(Base.read, (Seq.FASTAParser{Seq.AminoAcidSequence},))

precompile(Base.open, (ASCIIString, Type{Seq.FASTQ},))
precompile(Base.read, (Seq.FASTQParser{Seq.DNASequence},))

precompile(Base.open, (ASCIIString, Type{Intervals.BED},))
precompile(Base.read, (Intervals.BEDParser,))
