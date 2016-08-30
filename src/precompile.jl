# Precompile
# ==========

# Reader
# ------

precompile(Base.open, (String, Type{Seq.FASTA},))
precompile(Base.open, (String, Type{Seq.FASTQ},))
precompile(Base.open, (String, Type{Intervals.BED},))
precompile(Base.read, (Seq.FASTAReader{Seq.BioSequence},))
precompile(Base.read, (Seq.FASTAReader{Seq.DNASequence},))
precompile(Base.read, (Seq.FASTAReader{Seq.RNASequence},))
precompile(Base.read, (Seq.FASTAReader{Seq.AminoAcidSequence},))
precompile(Base.read, (Seq.FASTQReader{Seq.DNASequence},))
precompile(Base.read, (Intervals.BEDReader,))
