export getcoordarray,
    getrmsd,
    getdisps,
    getdist


import Base.-


"""Get the atomic coordinates of a `StrucElementOrList` as an array with each
column corresponding to one atom."""
function getcoordarray(element::StrucElementOrList, args...)
    atoms = collectatoms(element, args...)
    coords = zeros(3, length(atoms))
    for j in eachindex(atoms)
        coords[1,j] = getx(atoms[j])
        coords[2,j] = gety(atoms[j])
        coords[3,j] = getz(atoms[j])
    end
    return coords
end

getcoordarray(coords::Array{Float64}, args...) = coords


"""Get the root-mean-square deviation (RMSD) between two `StrucElementOrList`s
or coordinate arrays. Assumes they are aready aligned."""
function getrmsd(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate RMSD"
    diff = coords_one - coords_two
    return sqrt(sum(diff .* diff) / size(coords_one, 2))
end

# Repeat assertion here?
# Backbone by default?
getrmsd(element_one::StrucElementOrList, element_two::StrucElementOrList, args...) = getrmsd(getcoordarray(element_one, args...), getcoordarray(element_two, args...))


"""Get the displacements between atomic coordinates from two
`StrucElementOrList`s or coordinate arrays. Assumes they are aready aligned."""
function getdisps(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate displacements"
    diff = coords_one - coords_two
    return sqrt(sum(diff .* diff, 1))[:]
end

getdisps(element_one::StrucElementOrList, element_two::StrucElementOrList, args...) = getdisps(getcoordarray(element_one, args...), getcoordarray(element_two, args...))


# Maximum keyword to get max instead?
"""Get the minimum distance between two `StrucElementOrList`s or coordinate
arrays."""
function getdist(element_one::StrucElementOrList, element_two::StrucElementOrList, args...)
    coords_one = getcoordarray(element_one, args...)
    coords_two = getcoordarray(element_two, args...)
    count_one = size(coords_one, 2)
    count_two = size(coords_two, 2)
    sq_dists = zeros(count_one, count_two)
    for i in 1:count_one
        for j in 1:count_two
            @inbounds sq_dists[i,j] = (coords_one[1,i] - coords_two[1,j]) ^ 2 + (coords_one[2,i] - coords_two[2,j]) ^ 2 + (coords_one[3,i] - coords_two[3,j]) ^ 2
        end
    end
    return sqrt(minimum(sq_dists))
end

getdist(atom_one::AbstractAtom, atom_two::AbstractAtom) = sqrt((getx(atom_one) - getx(atom_two)) ^ 2 + (gety(atom_one) - gety(atom_two)) ^ 2 + (getz(atom_one) - getz(atom_two)) ^ 2)

-(element_one::StrucElementOrList, element_two::StrucElementOrList) = getdist(element_one, element_two)
