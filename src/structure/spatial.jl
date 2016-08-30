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
Additional arguments are atom selector functions - only atoms that return
`True` from the functions are retained.
"""
function coordarray(el::StructuralElementOrList, atom_selectors::Function...)
    atom_list = collectatoms(el, atom_selectors...)
    coords_out = zeros(3, length(atom_list))
    for j in eachindex(atom_list)
        coords_out[1,j] = x(atom_list[j])
        coords_out[2,j] = y(atom_list[j])
        coords_out[3,j] = z(atom_list[j])
    end
    return coords_out
end

# Selector functions ignored
coordarray(coords_in::Array{Float64}, atom_selectors::Function...) = coords_in


"""
Get the root-mean-square deviation (RMSD) between two `StructuralElementOrList`s
or two coordinate `Array`s of the same size. Assumes they are already
superimposed.
Additional arguments are atom selector functions - only atoms that return
`True` from the functions are retained.
"""
function rmsd(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate RMSD"
    diff = coords_one - coords_two
    return sqrt(sum(diff .* diff) / size(coords_one, 2))
end

rmsd(el_one::StructuralElementOrList, el_two::StructuralElementOrList, atom_selectors::Function...) = rmsd(coordarray(el_one, atom_selectors...), coordarray(el_two, atom_selectors...))


"""
Get the displacements between atomic coordinates from two
`StructuralElementOrList`s or two coordinate `Array`s of the same size. Assumes
they are already superimposed.
Additional arguments are atom selector functions - only atoms that return
`True` from the functions are retained.
"""
function displacements(coords_one::Array{Float64}, coords_two::Array{Float64})
    @assert size(coords_one) == size(coords_two) "Sizes of coordinate arrays are different - cannot calculate displacements"
    diff = coords_one - coords_two
    return sqrt(sum(diff .* diff, 1))[:]
end

displacements(el_one::StructuralElementOrList, el_two::StructuralElementOrList, atom_selectors::Function...) = displacements(coordarray(el_one, atom_selectors...), coordarray(el_two, atom_selectors...))


"""
Get the minimum square distance between two `StructuralElementOrList`s.
Additional arguments are atom selector functions - only atoms that return
`True` from the functions are retained.
"""
function sqdistance(el_one::StructuralElementOrList, el_two::StructuralElementOrList, atom_selectors::Function...)
    coords_one = coordarray(el_one, atom_selectors...)
    coords_two = coordarray(el_two, atom_selectors...)
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


"""
Get the minimum distance between two `StructuralElementOrList`s.
Additional arguments are atom selector functions - only atoms that return
`True` from the functions are retained.
"""
distance(el_one::StructuralElementOrList, el_two::StructuralElementOrList, atom_selectors::Function...) = sqrt(sqdistance(
        el_one,
        el_two,
        atom_selectors...
    ))

distance(atom_one::AbstractAtom, atom_two::AbstractAtom) = sqrt(sqdistance(atom_one, atom_two))


# in docs say we assume they are in a bond, and degrees/rads
bondangle(atom_a::AbstractAtom, atom_b::AbstractAtom, atom_c::AbstractAtom) = bondangle(
        coords(atom_a) - coords(atom_b),
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


"""

Additional arguments are residue selector functions - only residues that return
`True` from the functions are retained.
"""
function ramachandranangles(el::StructuralElementOrList, residue_selectors::Function...)
    res_list = collectresidues(el, residue_selectors...)
    @assert length(res_list) > 2 "Multiple residues required to calculate Ramachandran angles"
    phi_angles = Float64[]
    psi_angles = Float64[]
    for i in 2:length(res_list)-1
        res = res_list[i]
        res_prev = res_list[i-1]
        res_next = res_list[i+1]
        if sequentialresidues(res_prev, res) && sequentialresidues(res, res_next)
            # Eliminate this try
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


function contactmap(el_one::StructuralElementOrList, el_two::StructuralElementOrList, contact_dist::Real)
    sq_contact_dist = contact_dist ^ 2
    contacts = falses(length(el_one), length(el_two))
    el_one_list = collect(el_one)
    el_two_list = collect(el_two)
    for i in 1:length(el_one)
        for j in 1:length(el_two)
            if sqdistance(el_one_list[i], el_two_list[j]) <= sq_contact_dist
                contacts[i,j] = true
            end
        end
    end
    return contacts
end

function contactmap(el::StructuralElementOrList, contact_dist::Real)
    sq_contact_dist = contact_dist ^ 2
    contacts = eye(Bool, length(el))
    el_list = collect(el)
    for i in 1:length(el)
        for j in 1:i-1
            if sqdistance(el_list[i], el_list[j]) <= sq_contact_dist
                contacts[i,j] = true
                contacts[j,i] = true
            end
        end
    end
    return contacts
end
