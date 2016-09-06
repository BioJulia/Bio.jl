##BLAST Tools module for Bio.jl - Spec and Roadmap

Bio.jl needs a way to run BLAST and other command line tools commonly used in biology applications, taking inputs and capturing outputs in Bio.jl formats.

###Minimum requirements for BLAST module
With BLAST+ installed on users' system:

**Input**

- [x] accept seq object or file name as query and subject
    - [x] accept single sequence objects
    - [x] accept multi-sequence objects
- [x] blastn, blastp, ~~blastall~~
    - [ ] others?
- [x] accept flags for search modification
- [x] take/parse blast xml output

**Output**

- [x] return alignment object or strings for alignment
- [x] return stats object w/other information
- [ ] allow return of other blast outputs (tsv, csv etc)

###Other Useful Features

- [ ] bundle blast package from within Bio.jl (is this possible?)
- [ ] pull queries/subjects from NCBI
    - [ ] incorporate into Bio.jl formats
- [ ] run BLAST through NCBI servers instead of locally
- [ ] ...
