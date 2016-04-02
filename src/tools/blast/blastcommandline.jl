# wrapper for BLAST+ command line functions

immutable BlastResult
    bitscore::Float64
    expect::Float64
    queryname::ASCIIString
    hitname::ASCIIString
    hit::BioSequence
    alignment::AlignedSequence
end

"""
`readblastXML(blastrun::ASCIIString)`
Parse XML output of a blast run. Input is an XML string eg:
```julia
results = readall(open("blast_results.xml")) # need to use `readstring` instead of `readall` for v0.5
readblastXML(results)
```

Returns Vector{BlastResult} with the sequence of the hit, the Alignment with query sequence, bitscore and expect value
"""
function readblastXML(blastrun::ASCIIString; seqtype="nucl")
    xdoc = parse_string(blastrun)
    xroot = root(xdoc)
    params = get_elements_by_tagname(xroot, "BlastOutput_param")
    iterations = get_elements_by_tagname(xroot, "BlastOutput_iterations")
    results = BlastResult[]
    for iteration in collect(child_elements(iterations[1]))
        queryname = content(find_element(iteration, "Iteration_query-def"))
        hits = get_elements_by_tagname(iteration, "Iteration_hits")
        for hit in collect(child_elements(hits[1]))
            hitname = content(find_element(hit, "Hit_def"))
            hsps = get_elements_by_tagname(hit, "Hit_hsps")
            for hsp in collect(child_elements(hsps[1]))
                if seqtype == "nucl"
                    qseq = BioSequence{DNAAlphabet{4}}(content(find_element(hsp, "Hsp_qseq")))
                    hseq = BioSequence{DNAAlphabet{4}}(content(find_element(hsp, "Hsp_hseq")))
                elseif seqtype == "prot"
                    qseq = AminoAcidSequence(content(find_element(hsp, "Hsp_qseq")))
                    hseq = AminoAcidSequence(content(find_element(hsp, "Hsp_hseq")))
                else
                    throw(error("Please use \"nucl\" or \"prot\" for seqtype"))
                end

                aln = AlignedSequence(qseq, hseq)
                bitscore = float(content(find_element(hsp, "Hsp_bit-score")))
                expect = float(content(find_element(hsp, "Hsp_evalue")))
                push!(results, BlastResult(bitscore, expect, queryname, hitname, hseq, aln))
            end
        end
    end
    return results
end

"""
`readblastXML(blastrun::Cmd)`
Parse command line blast query with XML output. Input is the blast command line command.

Returns Vector{BlastResult} with the sequence of the hit, the Alignment with query sequence, bitscore and expect value
"""
function readblastXML(blastrun::Cmd; seqtype="nucl")
    # need to use `readstring` instead of `readall` for v0.5
    readblastXML(readall(blastrun), seqtype=seqtype)
end


"""
`blastn(query, subject, flags...)``
Runs blastn on `query` against `subject`.
    Subjects and queries may be file names (as strings), BioSequence{DNAAlphabet{4}} type or
    Array of BioSequence{DNAAlphabet{4}}.
    May include optional `flag`s such as `["-perc_identity", 95,]`. Do not use `-outfmt`.
"""
function blastn(query::AbstractString, subject::AbstractString, flags=[]; db::Bool=false)
    if db
        results = readblastXML(`blastn -query $query -db $subject $flags -outfmt 5`)
    else
        results = readblastXML(`blastn -query $query -subject $subject $flags -outfmt 5`)
    end
    return results
end

function blastn(query::BioSequence{DNAAlphabet{4}}, subject::BioSequence{DNAAlphabet{4}}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastn(querypath, subjectpath, flags)
end

function blastn{S <: BioSequence{DNAAlphabet{4}}}(query::BioSequence{DNAAlphabet{4}}, subject::Vector{S}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastn(querypath, subjectpath, flags)
end

function blastn(query::BioSequence{DNAAlphabet{4}}, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        blastn(querypath, subject, flags, db=true)
    else
        blastn(querypath, subject, flags)
    end
end

function blastn{S <: BioSequence{DNAAlphabet{4}}}(query::Vector{S}, subject::Vector{S}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastn(querypath, subjectpath, flags)
end

"""
`blastp(query, subject, flags...)``
Runs blastn on `query` against `subject`.
    Subjects and queries may be file names (as strings), `BioSequence{AminoAcidSequence}` type or
    Array of `BioSequence{AminoAcidSequence}`.
    May include optional `flag`s such as `["-perc_identity", 95,]`. Do not use `-outfmt`.
"""
function blastp(query::AbstractString, subject::AbstractString, flags=[]; db::Bool=false)
    if db
        results = readblastXML(`blastp -query $query -db $subject $flags -outfmt 5`, seqtype = "prot")
    else
        results = readblastXML(`blastp -query $query -subject $subject $flags -outfmt 5`, seqtype = "prot")
    end
end

function blastp(query::AminoAcidSequence, subject::AminoAcidSequence, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastp(querypath, subjectpath, flags)
end

function blastp{S <: AminoAcidSequence}(query::AminoAcidSequence, subject::Vector{S}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastp(querypath, subjectpath, flags)
end

function blastp(query::AminoAcidSequence, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        blastp(querypath, subject, flags, db=true)
    else
        blastp(querypath, subject, flags)
    end
end

function blastp{S <: AminoAcidSequence}(query::Vector{S}, subject::Vector{S}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastp(querypath, subjectpath, flags)
end

"""
`makefasta(sequence::BioSequence)`
Create temporary fasta-formated file for blasting.
"""
function makefasta{T <: BioSequence}(sequence::T)
    path, io = mktemp()
    write(io, ">$path\n$(convert(AbstractString, sequence))\n")
    close(io)
    return path
end

"""
`makefasta(sequence::Vector{BioSequence})`
Create temporary multi fasta-formated file for blasting.
"""
function makefasta{T <: BioSequence}(sequences::Vector{T})
    path, io = mktemp()
    counter = 1
    for sequence in sequences
        write(io, ">$path$counter\n$(convert(AbstractString, sequence))\n")
        counter += 1
    end
    close(io)
    return path
end
