export coordarray,
    rmsd,
    disps,
    dist


import Base.-


"""Get the atomic coordinates of a `StrucElementOrList` as an array with each
column corresponding to one atom."""
function coordarray(element::StrucElementOrList, selector_functions::Function...)
    atoms = collectatoms(element, selector_functions...)
    coords = zeros(3, length(atoms))
    for j in eachindex(atoms)
        coords[1,j] = x(atoms[j])
        coords[2,j] = y(atoms[j])
        coords[3,j] = z(atoms[j])
    end
    return coords
end

coordarray(coords::Array{Float64}, selector_functions::Function...) = coords


"""Get the root-mean-square deviation (RMSD) between two `StrucElementOrList`s
or coordinate arrays. Assumes they are aready aligned."""
function rmsd(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate RMSD"
    diff = coords_one - coords_two
    return sqrt(sum(diff .* diff) / size(coords_one, 2))
end

# Repeat assertion here?
# Backbone by default?
rmsd(element_one::StrucElementOrList, element_two::StrucElementOrList, selector_functions::Function...) = rmsd(coordarray(element_one, selector_functions...), coordarray(element_two, selector_functions...))


"""Get the displacements between atomic coordinates from two
`StrucElementOrList`s or coordinate arrays. Assumes they are aready aligned."""
function disps(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate displacements"
    diff = coords_one - coords_two
    return sqrt(sum(diff .* diff, 1))[:]
end

disps(element_one::StrucElementOrList, element_two::StrucElementOrList, selector_functions::Function...) = disps(coordarray(element_one, selector_functions...), coordarray(element_two, selector_functions...))


"""Get the minimum distance between two `StrucElementOrList`s or coordinate
arrays."""
function dist(element_one::StrucElementOrList, element_two::StrucElementOrList, selector_functions::Function...)
    coords_one = coordarray(element_one, selector_functions...)
    coords_two = coordarray(element_two, selector_functions...)
    count_one = size(coords_one, 2)
    count_two = size(coords_two, 2)
    sq_dists = zeros(count_one, count_two)
    for i in 1:count_one, j in 1:count_two
        @inbounds sq_dists[i,j] = (coords_one[1,i] - coords_two[1,j]) ^ 2 + (coords_one[2,i] - coords_two[2,j]) ^ 2 + (coords_one[3,i] - coords_two[3,j]) ^ 2
    end
    return sqrt(minimum(sq_dists))
end

dist(atom_one::AbstractAtom, atom_two::AbstractAtom) = sqrt((x(atom_one) - x(atom_two)) ^ 2 + (y(atom_one) - y(atom_two)) ^ 2 + (z(atom_one) - z(atom_two)) ^ 2)

-(element_one::StrucElementOrList, element_two::StrucElementOrList) = dist(element_one, element_two)
