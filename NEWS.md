Bio.jl v0.5.0 Release Notes
===========================

* `Nucleotide`, `DNANucleotide` and `RNANucleotide` are renamed to `NucleicAcid`, `DNA` and `RNA`, respectively ([#391]).
* `qualities(record::SAMRecord)` returns decoded base qualities as `Vector{UInt8}` instead of `String` by default; use `qualities(String, record::SAMRecord)` to get a `String` value ([#399]).

[#391]: https://github.com/BioJulia/Bio.jl/issues/391
[#399]: https://github.com/BioJulia/Bio.jl/pull/399
