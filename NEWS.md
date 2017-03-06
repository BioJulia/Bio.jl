Bio.jl v0.5.0 Release Notes
===========================

* `Nucleotide`, `DNANucleotide` and `RNANucleotide` are renamed to `NucleicAcid`, `DNA` and `RNA`, respectively ([#391]).
* The SAM reader is redesigned from scratch to improve the performance ([#399]).
* `sequence(record::SAMRecord)` returns `DNASequence` instead of `String` by default; use `sequence(String, record)` to get a `String` value ([#399]).
* `qualities(record::SAMRecord)` returns decoded base qualities as `Vector{UInt8}` instead of `String` by default; use `qualities(String, record)` to get a `String` value ([#399]).
* Support VCF and BCF file formats ([#378]).
* Support APIs for Entrez Programming Utilities ([#350]).
* Support GFF3 file format ([#138]).

[#138]: https://github.com/BioJulia/Bio.jl/pull/138
[#350]: https://github.com/BioJulia/Bio.jl/pull/350
[#378]: https://github.com/BioJulia/Bio.jl/pull/378
[#391]: https://github.com/BioJulia/Bio.jl/issues/391
[#399]: https://github.com/BioJulia/Bio.jl/pull/399
