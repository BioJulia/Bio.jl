# FASTA sequence types

immutable FASTA <: FileFormat end


"Metadata for FASTA sequence records containing just a `description` field"
type FASTAMetadata
    description::StringField
end


function FASTAMetadata()
    return FASTAMetadata(StringField())
end


function Base.(:(==))(a::FASTAMetadata, b::FASTAMetadata)
    return a.description == b.description
end


function Base.copy(metadata::FASTAMetadata)
    return FASTAMetadata(copy(metadata.description))
end


"FASTASeqRecord{S} is a `SeqRecord` for FASTA sequences of type `S`"
typealias FASTASeqRecord          SeqRecord{Sequence, FASTAMetadata}

"A `SeqRecord` type for FASTA DNA sequences"
typealias FASTADNASeqRecord       DNASeqRecord{FASTAMetadata}

"A `SeqRecord` type for FASTA RNA sequences"
typealias FASTARNASeqRecord       RNASeqRecord{FASTAMetadata}

"A `SeqRecord` type for FASTA amino acid sequences"
typealias FASTAAminoAcidSeqRecord AminoAcidSeqRecord{FASTAMetadata}

function Base.show{S}(io::IO, seqrec::SeqRecord{S, FASTAMetadata})
    write(io, ">", seqrec.name, " ", seqrec.metadata.description, "\n")
    show(io, seqrec.seq)
end

function Base.print{S}(io::IO, seqrec::SeqRecord{S,FASTAMetadata})
    write(io, ">", seqrec.name, " ", seqrec.metadata.description, "\n")
    print(io, seqrec.seq)
end


"Writes a FASTASeqRecord to an IO-stream (and obeys FASTAs max character constraint)"
function Base.write{T}(io::IO, seqrec::SeqRecord{T, FASTAMetadata})
    write(io, ">", seqrec.name)
    if !isempty(seqrec.metadata.description)
        write(io, " ", seqrec.metadata.description)
    end
    write(io, "\n")
    maxchars = 79
    counter = 1
    len = length(seqrec.seq)
    for nt in seqrec.seq
        show(io, nt)
        if counter % maxchars == 0 && counter < len
            write(io, "\n")
        end
        counter += 1
    end
    write(io, "\n")
end
