

# Genotype Array, a type designed to store genotype counts.
# Behaves like a named array, and wraps a PooledDataArray.
# Allows subsetting by Loci names, and by individual names.

# Data is a sparse array containing allele counts.
# AllelesPerLocus tells how many alleles are at each locus.
# LocusFactor describes which columns of the matrix, are for each locus.
# LocusAlleleNames has the names of the alleles for each locus.

type GenotypeArray{S <: AbstractString, N}
    Individuals::Vector{S}
    LocusFactor::PooledDataArray{S}
    LocusAlleleNames::Dict{S, Vector{S}}
    Data::PooledDataArray{Int, N}
    #AllelesPerLocus::Vector{Int}
    #Ploidy::Vector{Int}
end

# GenotypeVector has only one row - one individual.
typealias GenotypeVector{S} GenotypeArray{S, 1}
# GenotypeMatrix has multiple rows - multiple individuals.
typealias GenotypeMatrix{S} GenotypeArray{S, 2}









#=

IGNORE FOR NOW - SOME NON FUNCTIONAL CODE TO DO WITH
STORING GENOTYPED SNPS IN BINARY FORMAT. B.W.


immutable LocusSpec
    major::Char
    minor::Char
end


type Genotypes{C, P}
    genomes::Vector{BinarySNPs}
    numberOfLoci::Int
    individualNames::Vector{ASCIIString}
    lociNames::Vector{ASCIIString}
    loci::Vector{LocusSpec}
    chromosome::PooledDataArray{C}
    position::Vector{Int}
    ploidy::Vector{Int}
    population::PooledDataArray{P}
    strata::Nullable{DataFrame}
    hierarchy# Nullable Forumula - find formula class used by JuliaStats.
end


function Genotypes{C, P}(hier, snpbins::Vector{BinarySNPs} = BinarySNPs[],
                         indnames::Vector{ASCIIString} = ASCIIString[],
                         locinames::Vector{ASCIIString} = ASCIIString[],
                         loc::Vector{LocusSpec} = LocusSpec[],
                         chrom::PooledDataArray{C} = PooledDataArray(C),
                         pos::Vector{Int} = Int[],
                         pop::PooledDataArray{P} = PooledDataArray(P),
                         strat::DataFrame = DataFrame(),
                         ploid::Vector{Int} = Int[])

    if length(genomes) > 0
        # Check the number of loci in each genome provided, and the one provided.
        nloci = genomes[1].numberOfLoci
        for snpbin in snpbins
            if snpbin.numberOfLoci != nloc
                error("The number of loci in each input genome is not equal.")
            end
        end
        # Check the labels provided for individuals.
        if length(indnames) == 0
            for snpbin in snpbins
                push!(indnames, snpbin.label)
            end
        elseif length(indnames) != length(genomes)
            error("The number of individual names, and genomes provided are inconsistent.")
        end
        # Check ploidy provided.
        if length(ploidy) == 0
            for snpbin in snpbins
                push!(ploidy, genome.ploidy)
            end
        end
        # Check loci names provided.
        if length(locinames) > 0 && length(locinames) != nloc
            error("The number of loci names provided is inconsistent with the number of loci.")
        end
        # Check loc provided.
        if length(loc) > 0 && length(loc) != nloc
            error("The number of loci provided is inconsistent with the number of loci.")
        end
        # Check the chrom argument.
        if length(chrom) > 0 && length(chrom) != nloc
            error("The number of chromosome annotations is inconsistent with the number of loci.")
        end
        # Check length of the position argument.
        if length(pos) > 0 && length(pos) != nloc
            error("The number of positions is inconsistent with the number of loci.")
        end
        # Check the length of the pop argument.
        if length(pop) > 0 && length(pop) != length(genomes)
            error("The number of population annotations is inconsistent with the number of genomes.")
        end
        # Check the strata argument.



    end

    return Genotypes{C, P}(genomes, nloc, indnames, locnames, loc, chrom, pos, ploid, pop,
               strat, hier)
end

=#
