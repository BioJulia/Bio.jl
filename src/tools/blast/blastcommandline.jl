# BLAST+ Wrapper
# ==============
#
# Wrapper for BLAST+ command line functions.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable BLASTResult
    bitscore::Float64
    expect::Float64
    queryname::String
    hitname::String
    hit::BioSequence
    alignment::AlignedSequence
end

"""
    readblastXML(blastrun::AbstractString)

Parse XML output of a blast run. Input is an XML string eg:
```julia
results = readstring(open("blast_results.xml"))
readblastXML(results)
```

Returns Vector{BLASTResult} with the sequence of the hit, the Alignment with query sequence, bitscore and expect value
"""
function readblastXML(blastrun::AbstractString; seqtype="nucl")
    dc = parsexml(blastrun)
    rt = root(dc)
    results = BLASTResult[]
    for iteration in find(rt, "/BlastOutput/BlastOutput_iterations/Iteration")
        queryname = content(findfirst(iteration, "Iteration_query-def"))
        for hit in find(iteration, "Iteration_hits")
            if countelements(hit) > 0
                hitname = content(findfirst(hit, "./Hit/Hit_def"))
                hsps = findfirst(hit, "./Hit/Hit_hsps")
                if seqtype == "nucl"
                    qseq = DNASequence(content(findfirst(hsps, "./Hsp/Hsp_qseq")))
                    hseq = DNASequence(content(findfirst(hsps, "./Hsp/Hsp_hseq")))
                elseif seqtype == "prot"
                    qseq = AminoAcidSequence(content(findfirst(hsps, "./Hsp/Hsp_qseq")))
                    hseq = AminoAcidSequence(content(findfirst(hsps, "./Hsp/Hsp_hseq")))
                else
                    throw(error("Please use \"nucl\" or \"prot\" for seqtype"))
                end

                aln = AlignedSequence(qseq, hseq)
                bitscore = float(content(findfirst(hsps, "./Hsp/Hsp_bit-score")))
                expect = float(content(findfirst(hsps, "./Hsp/Hsp_evalue")))
                push!(results, BLASTResult(bitscore, expect, queryname, hitname, hseq, aln))
            end
        end
    end
    return results
end

"""
`readblastXML(blastrun::Cmd)`
Parse command line blast query with XML output. Input is the blast command line command, eg:
```julia
blastresults = `blastn -query seq1.fasta -db some_database -outfmt 5`
readblastXML(blastresults)
```

Returns Vector{BLASTResult} with the sequence of the hit, the Alignment with query sequence, bitscore and expect value
"""
function readblastXML(blastrun::Cmd; seqtype="nucl")
    return readblastXML(readstring(blastrun), seqtype=seqtype)
end


"""
`blastn(query, subject, flags...)``
Runs blastn on `query` against `subject`.
    Subjects and queries may be file names (as strings), DNASequence type or
    Array of DNASequence.
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

function blastn(query::DNASequence, subject::DNASequence, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastn(querypath, subjectpath, flags)
end

function blastn(query::DNASequence, subject::Vector{DNASequence}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    blastn(querypath, subjectpath, flags)
end

function blastn(query::DNASequence, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        return blastn(querypath, subject, flags, db=true)
    else
        return blastn(querypath, subject, flags)
    end
end

function blastn(query::Vector{DNASequence}, subject::Vector{DNASequence}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastn(querypath, subjectpath, flags)
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
    return results
end

function blastp(query::AminoAcidSequence, subject::AminoAcidSequence, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastp(querypath, subjectpath, flags)
end

function blastp(query::AminoAcidSequence, subject::Vector{AminoAcidSequence}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastp(querypath, subjectpath, flags)
end

function blastp(query::AminoAcidSequence, subject::AbstractString, flags=[]; db::Bool=false)
    querypath = makefasta(query)
    if db
        return blastp(querypath, subject, flags, db=true)
    else
        return blastp(querypath, subject, flags)
    end
end

function blastp(query::Vector{AminoAcidSequence}, subject::Vector{AminoAcidSequence}, flags=[])
    querypath, subjectpath = makefasta(query), makefasta(subject)
    return blastp(querypath, subjectpath, flags)
end

# Create temporary fasta-formated file for blasting.
function makefasta(sequence::BioSequence)
    path, io = mktemp()
    write(io, ">$path\n$(convert(String, sequence))\n")
    close(io)
    return path
end

# Create temporary multi fasta-formated file for blasting.
function makefasta{T <: BioSequence}(sequences::Vector{T})
    path, io = mktemp()
    counter = 1
    for sequence in sequences
        write(io, ">$path$counter\n$(convert(String, sequence))\n")
        counter += 1
    end
    close(io)
    return path
end
