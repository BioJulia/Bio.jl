export getcoordarray,
    getrmsd,
    getdisps,
    getdist


import Base.-


# This doesn't work yet
getcoordarray(element::StrucElementOrList, args...) = map(getcoords, collectatoms(element, args...))


function getrmsd(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate RMSD"
    diff = coords_one - coords_two
    return sqrt(sum(diff .* diff) / size(coords_one, 2))
end

# Repeat assertion here?
getrmsd(element_one::StrucElementOrList, element_two::StrucElementOrList, args...) = getrmsd(getcoordarray(element_one, args...), getcoordarray(element_two, args...))


function getdisps(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate displacements"
    diff = coords_one - coords_two
    return sqrt(sum(diff .* diff, 1))[:]
end

getdisps(element_one::StrucElementOrList, element_two::StrucElementOrList, args...) = getdisps(getcoordarray(element_one, args...), getcoordarray(element_two, args...))


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
