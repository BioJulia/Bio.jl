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
`true` from the functions are retained.
"""
function coordarray(el::StructuralElementOrList, atom_selectors::Function...)
    at_list = collectatoms(el, atom_selectors...)
    coords_out = zeros(3, length(at_list))
    for j in eachindex(at_list)
        coords_out[1,j] = x(at_list[j])
        coords_out[2,j] = y(at_list[j])
        coords_out[3,j] = z(at_list[j])
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
`true` from the functions are retained.
"""
function rmsd(coords_one::Array{Float64}, coords_two::Array{Float64})
    if size(coords_one) != size(coords_two)
        throw(ArgumentError("Sizes of coordinate arrays are different - cannot calculate RMSD"))
    end
    diff = coords_one - coords_two
    return sqrt.(dot(diff, diff) / size(coords_one, 2))
end

function rmsd(el_one::StructuralElementOrList,
            el_two::StructuralElementOrList,
            atom_selectors::Function...)
    return rmsd(
        coordarray(el_one, atom_selectors...),
        coordarray(el_two, atom_selectors...))
end


"""
Get the displacements between atomic coordinates from two
`StructuralElementOrList`s or two coordinate `Array`s of the same size. Assumes
they are already superimposed.
Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function displacements(coords_one::Array{Float64}, coords_two::Array{Float64})
    if size(coords_one) != size(coords_two)
        throw(ArgumentError("Sizes of coordinate arrays are different - cannot calculate displacements"))
    end
    diff = coords_one - coords_two
    return sqrt.(sum(diff .* diff, 1))[:]
end

function displacements(el_one::StructuralElementOrList,
                    el_two::StructuralElementOrList,
                    atom_selectors::Function...)
    return displacements(
        coordarray(el_one, atom_selectors...),
        coordarray(el_two, atom_selectors...))
end


"""
Get the minimum square distance between two `StructuralElementOrList`s.
Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function sqdistance(el_one::StructuralElementOrList,
                    el_two::StructuralElementOrList,
                    atom_selectors::Function...)
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

function sqdistance(at_one::AbstractAtom, at_two::AbstractAtom)
    return (x(at_one) - x(at_two)) ^ 2 +
           (y(at_one) - y(at_two)) ^ 2 +
           (z(at_one) - z(at_two)) ^ 2
end


"""
Get the minimum distance between two `StructuralElementOrList`s.
Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function Bio.distance(el_one::StructuralElementOrList,
                      el_two::StructuralElementOrList,
                      atom_selectors::Function...)
    return sqrt(sqdistance(el_one, el_two, atom_selectors...))
end

function Bio.distance(at_one::AbstractAtom, at_two::AbstractAtom)
    return sqrt(sqdistance(at_one, at_two))
end


"""
Calculate the bond angle in radians between three `AbstractAtom`s, where A-B and
B-C are assumed to be bonded.
"""
function bondangle(at_a::AbstractAtom,
                at_b::AbstractAtom,
                at_c::AbstractAtom)
    return bondangle(
        coords(at_a) - coords(at_b),
        coords(at_c) - coords(at_b)
    )
end

function bondangle(vec_a::Vector{Float64}, vec_b::Vector{Float64})
    return acos(dot(vec_a, vec_b) / (norm(vec_a) * norm(vec_b)))
end


"""
Calculate the dihedral angle in radians defined by four `AbstractAtom`s. This is
the angle between the planes defined by atoms (A,B,C) and (B,C,D).
"""
function dihedralangle(at_a::AbstractAtom,
            at_b::AbstractAtom,
            at_c::AbstractAtom,
            at_d::AbstractAtom)
    return dihedralangle(
        coords(at_b) - coords(at_a),
        coords(at_c) - coords(at_b),
        coords(at_d) - coords(at_c))
end

function dihedralangle(vec_a::Vector{Float64},
                    vec_b::Vector{Float64},
                    vec_c::Vector{Float64})
    return atan2(
        dot(cross(cross(vec_a, vec_b), cross(vec_b, vec_c)), vec_b / norm(vec_b)),
        dot(cross(vec_a, vec_b), cross(vec_b, vec_c)))
end


"""
Calculate the omega angle in radians for an `AbstractResidue`. Arguments are the
residue (atoms "N" and "CA" required) and the previous residue (atoms "CA" and
"C" required).
"""
function omegaangle(res::AbstractResidue, res_prev::AbstractResidue)
    if !("CA" in atomnames(res_prev))
        throw(ArgumentError("Atom with atom name \"CA\" not found in previous residue"))
    elseif !("C"  in atomnames(res_prev))
        throw(ArgumentError("Atom with atom name \"C\" not found in previous residue"))
    elseif !("N"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in atomnames(res))
        throw(ArgumentError("Atom with atom name \"CA\" not found in residue"))
    end
    return dihedralangle(res_prev["CA"], res_prev["C"], res["N"], res["CA"])
end

"""
Calculate the phi angle in radians for an `AbstractResidue`. Arguments are the
residue (atoms "N", "CA" and "C" required) and the previous residue (atom "C"
required).
"""
function phiangle(res::AbstractResidue, res_prev::AbstractResidue)
    if !("C"  in atomnames(res_prev))
        throw(ArgumentError("Atom with atom name \"C\" not found in previous residue"))
    elseif !("N"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in atomnames(res))
        throw(ArgumentError("Atom with atom name \"CA\" not found in residue"))
    elseif !("C"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"C\" not found in residue"))
    end
    return dihedralangle(res_prev["C"], res["N"], res["CA"], res["C"])
end

"""
Calculate the psi angle in radians for an `AbstractResidue`. Arguments are the
residue (atoms "N", "CA" and "C" required) and the next residue (atom "N"
required).
"""
function psiangle(res::AbstractResidue, res_next::AbstractResidue)
    if !("N"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in atomnames(res))
        throw(ArgumentError("Atom with atom name \"CA\" not found in residue"))
    elseif !("C"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"C\" not found in residue"))
    elseif !("N"  in atomnames(res_next))
        throw(ArgumentError("Atom with atom name \"N\" not found in next residue"))
    end
    return dihedralangle(res["N"], res["CA"], res["C"], res_next["N"])
end


"""
Calculate the `Vector`s of phi and psi angles of a `StructuralElementOrList`.
The vectors have `NaN` for residues where an angle cannot be calculated,
e.g. due to missing atoms or lack of an adjacent residue.
Additional arguments are residue selector functions - only residues that return
`true` from the functions are retained.
"""
function ramachandranangles(el::StructuralElementOrList,
                    residue_selectors::Function...)
    res_list = collectresidues(el, residue_selectors...)
    if length(res_list) < 2
        throw(ArgumentError("Multiple residues required to calculate Ramachandran angles"))
    end
    phi_angles = Float64[NaN] # First res has no previous res
    psi_angles = Float64[]
    # Phi angles
    for i in 2:length(res_list)
        res = res_list[i]
        res_prev = res_list[i-1]
        if sequentialresidues(res_prev, res)
            try
                phi_angle = phiangle(res, res_prev)
                push!(phi_angles, phi_angle)
            catch ex
                isa(ex, ArgumentError) || rethrow()
                push!(phi_angles, NaN)
            end
        else
            push!(phi_angles, NaN)
        end
    end
    # Psi angles
    for i in 1:length(res_list)-1
        res = res_list[i]
        res_next = res_list[i+1]
        if sequentialresidues(res, res_next)
            try
                psi_angle = psiangle(res, res_next)
                push!(psi_angles, psi_angle)
            catch ex
                isa(ex, ArgumentError) || rethrow()
                push!(psi_angles, NaN)
            end
        else
            push!(psi_angles, NaN)
        end
    end
    push!(psi_angles, NaN) # Last res has no next res
    return phi_angles, psi_angles
end


"""
Calculate the contact map for a `StructuralElementOrList`, or between two
`StructuralElementOrList`s. This is a `BitArray{2}` with `true` where the
sub-elements are no further than the contact distance and `false` otherwise.
"""
function contactmap(el_one::StructuralElementOrList,
                el_two::StructuralElementOrList,
                contact_dist::Real)
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
    contacts = falses(length(el), length(el))
    el_list = collect(el)
    for i in 1:length(el)
        contacts[i,i] = true
        for j in 1:i-1
            if sqdistance(el_list[i], el_list[j]) <= sq_contact_dist
                contacts[i,j] = true
                contacts[j,i] = true
            end
        end
    end
    return contacts
end
