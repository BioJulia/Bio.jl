export
    coordarray,
    rmsd,
    displacements,
    sqdistance,
    distance,
    bondangle,
    dihedralangle,
    omegaangle,
    phiangle,
    psiangle,
    ramachandranangles,
    contactmap


"""
Get the atomic coordinates of a `StructuralElementOrList` as a 2D `Array` with
each column corresponding to one atom.
"""
function coordarray(element::StructuralElementOrList, selector_functions::Function...)
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


"""
Get the root-mean-square deviation (RMSD) between two `StructuralElementOrList`s
or two coordinate `Array`s of the same size. Assumes they are already
superimposed.
"""
function rmsd(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate RMSD"
    difference = coords_one - coords_two
    return sqrt(sum(difference .* difference) / size(coords_one, 2))
end

rmsd(element_one::StructuralElementOrList, element_two::StructuralElementOrList, selector_functions::Function...) = rmsd(coordarray(element_one, selector_functions...), coordarray(element_two, selector_functions...))


"""
Get the displacements between atomic coordinates from two
`StructuralElementOrList`s or two coordinate `Array`s of the same size. Assumes
they are already superimposed.
"""
function displacements(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate displacements"
    difference = coords_one - coords_two
    return sqrt(sum(difference .* difference, 1))[:]
end

displacements(element_one::StructuralElementOrList, element_two::StructuralElementOrList, selector_functions::Function...) = displacements(coordarray(element_one, selector_functions...), coordarray(element_two, selector_functions...))


function sqdistance(element_one::StructuralElementOrList, element_two::StructuralElementOrList, selector_functions::Function...)
    coords_one = coordarray(element_one, selector_functions...)
    coords_two = coordarray(element_two, selector_functions...)
    min_sq_dist = Inf
    for i in 1:size(coords_one, 2)
        for j in 1:size(coords_two, 2)
            @inbounds sq_dist = (coords_one[1,i] - coords_two[1,j]) ^ 2 + (coords_one[2,i] - coords_two[2,j]) ^ 2 + (coords_one[3,i] - coords_two[3,j]) ^ 2
            if sq_dist < min_sq_dist
                min_sq_dist = sq_dist
            end
        end
    end
    return min_sq_dist
end

sqdistance(atom_one::AbstractAtom, atom_two::AbstractAtom) = (x(atom_one) - x(atom_two)) ^ 2 + (y(atom_one) - y(atom_two)) ^ 2 + (z(atom_one) - z(atom_two)) ^ 2


"Get the minimum distance between two `StructuralElementOrList`s."
distance(element_one::StructuralElementOrList, element_two::StructuralElementOrList, selector_functions::Function...) = sqrt(sqdistance(
        element_one,
        element_two,
        selector_functions...
    ))

distance(atom_one::AbstractAtom, atom_two::AbstractAtom) = sqrt(sqdistance(atom_one, atom_two))


# in docs say we assume they are in a bond, and degrees/rads
bondangle(atom_a::AbstractAtom, atom_b::AbstractAtom, atom_c::AbstractAtom) = bondangle(
        coords(atom_b) - coords(atom_a),
        coords(atom_c) - coords(atom_b)
    )
bondangle(vec_a::Vector{Float64}, vec_b::Vector{Float64}) = acos(
    dot(vec_a, vec_b) / (norm(vec_a) * norm(vec_b)))


# say that b and c are in common - angle between planes defined by (a,b,c) and (b,c,d)
dihedralangle(atom_a::AbstractAtom, atom_b::AbstractAtom, atom_c::AbstractAtom, atom_d::AbstractAtom) = dihedralangle(
        coords(atom_b) - coords(atom_a),
        coords(atom_c) - coords(atom_b),
        coords(atom_d) - coords(atom_c)
    )
dihedralangle(vec_a::Vector{Float64}, vec_b::Vector{Float64}, vec_c::Vector{Float64}) = atan2(
        dot(cross(cross(vec_a, vec_b), cross(vec_b, vec_c)), vec_b / norm(vec_b)), dot(cross(vec_a, vec_b), cross(vec_b, vec_c))
    )


# Be clear in docs that the argument order is not uniform
function omegaangle(res::AbstractResidue, res_prev::AbstractResidue)
    @assert haskey(atoms(res_prev), "CA") "Atom with atom name \"CA\" not found in previous residue"
    @assert haskey(atoms(res_prev), "C") "Atom with atom name \"C\" not found in previous residue"
    @assert haskey(atoms(res), "N") "Atom with atom name \"N\" not found in residue"
    @assert haskey(atoms(res), "CA") "Atom with atom name \"CA\" not found in residue"
    return dihedralangle(res_prev["CA"], res_prev["C"], res["N"], res["CA"])
end

function phiangle(res::AbstractResidue, res_prev::AbstractResidue)
    @assert haskey(atoms(res_prev), "C") "Atom with atom name \"C\" not found in previous residue"
    @assert haskey(atoms(res), "N") "Atom with atom name \"N\" not found in residue"
    @assert haskey(atoms(res), "CA") "Atom with atom name \"CA\" not found in residue"
    @assert haskey(atoms(res), "C") "Atom with atom name \"C\" not found in residue"
    return dihedralangle(res_prev["C"], res["N"], res["CA"], res["C"])
end

function psiangle(res::AbstractResidue, res_next::AbstractResidue)
    @assert haskey(atoms(res), "N") "Atom with atom name \"N\" not found in residue"
    @assert haskey(atoms(res), "CA") "Atom with atom name \"CA\" not found in residue"
    @assert haskey(atoms(res), "C") "Atom with atom name \"C\" not found in residue"
    @assert haskey(atoms(res_next), "N") "Atom with atom name \"N\" not found in next residue"
    return dihedralangle(res["N"], res["CA"], res["C"], res_next["N"])
end


function ramachandranangles(element::StructuralElementOrList, selector_functions::Function...)
    residues = collectresidues(element, selector_functions...)
    @assert length(residues) > 2 "Multiple residues required to calculate Ramachandran angles"
    phi_angles = Float64[]
    psi_angles = Float64[]
    for i in 2:length(residues)-1
        res = residues[i]
        res_prev = residues[i-1]
        res_next = residues[i+1]
        if sequentialresidues(res_prev, res) && sequentialresidues(res, res_next)
            try
                phi_angle = phiangle(res, res_prev)
                psi_angle = psiangle(res, res_next)
                push!(phi_angles, phi_angle)
                push!(psi_angles, psi_angle)
            end
        end
    end
    return phi_angles, psi_angles
end


function contactmap(element_one::StructuralElementOrList, element_two::StructuralElementOrList, contact_dist::Real)
    sq_contact_dist = contact_dist ^ 2
    contacts = falses(length(element_one), length(element_two))
    element_one_list = collect(element_one)
    element_two_list = collect(element_two)
    for i in 1:length(element_one)
        for j in 1:length(element_two)
            if sqdistance(element_one_list[i], element_two_list[j]) <= sq_contact_dist
                contacts[i,j] = true
            end
        end
    end
    return contacts
end

function contactmap(element::StructuralElementOrList, contact_dist::Real)
    sq_contact_dist = contact_dist ^ 2
    contacts = eye(Bool, length(element))
    element_list = collect(element)
    for i in 1:length(element)
        for j in 1:i-1
            if sqdistance(element_list[i], element_list[j]) <= sq_contact_dist
                contacts[i,j] = true
                contacts[j,i] = true
            end
        end
    end
    return contacts
end
