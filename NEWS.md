Bio.jl v0.5.0 Release Notes
===========================

* `Nucleotide`, `DNANucleotide` and `RNANucleotide` are renamed to `NucleicAcid`, `DNA` and `RNA`, respectively ([#391]).
* The SAM reader is redesigned from scratch to improve the performance ([#399]).
* `sequence(record::SAMRecord)` returns `DNASequence` instead of `String` by default; use `sequence(String, record)` to get a `String` value ([#399]).
* `qualities(record::SAMRecord)` returns decoded base qualities as `Vector{UInt8}` instead of `String` by default; use `qualities(String, record)` to get a `String` value ([#399]).
* Support VCF and BCF file formats ([#378]).
* Support APIs for Entrez Programming Utilities ([#350]).
* Support GFF3 file format ([#138]).
* Sequence names are ordered lexicographically by default ([#291]).
* IO APIs are reorganized into submodules:
    * `Bio.Seq.FASTAReader` => `Bio.Seq.FASTA.Reader` ([#413]).

[#138]: https://github.com/BioJulia/Bio.jl/pull/138
[#291]: https://github.com/BioJulia/Bio.jl/issues/291
[#350]: https://github.com/BioJulia/Bio.jl/pull/350
[#378]: https://github.com/BioJulia/Bio.jl/pull/378
[#391]: https://github.com/BioJulia/Bio.jl/issues/391
[#399]: https://github.com/BioJulia/Bio.jl/pull/399
[#413]: https://github.com/BioJulia/Bio.jl/pull/413
