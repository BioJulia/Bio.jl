export
    StructuralElement,
    AbstractAtom,
    Atom,
    DisorderedAtom,
    AbstractResidue,
    Residue,
    DisorderedResidue,
    Chain,
    Model,
    ProteinStructure,
    StructuralElementOrList,
    ishetatom,
    serial,
    atomname,
    altlocid,
    resname,
    chainid,
    resnumber,
    inscode,
    x,
    y,
    z,
    coords,
    occupancy,
    tempfac,
    element,
    charge,
    ishetero,
    isdisorderedatom,
    x!,
    y!,
    z!,
    coords!,
    defaultaltlocid,
    defaultatom,
    altlocids,
    resid,
    ishetres,
    atomnames,
    atoms,
    isdisorderedres,
    disorderedres,
    defaultresname,
    defaultresidue,
    resnames,
    resids,
    residues,
    modelnumber,
    chainids,
    chains,
    structurename,
    modelnumbers,
    models,
    defaultmodel,
    sequentialresidues,
    applyselectors,
    applyselectors!,
    collectresidues,
    collectatoms,
    countmodels,
    countchains,
    countresidues,
    countatoms,
    organise,
    formatomlist,
    choosedefaultaltlocid,
    organisemodel,
    organisestructure,
    stdatomselector,
    hetatomselector,
    atomnameselector,
    calpha_atom_names,
    calphaselector,
    cbeta_atom_names,
    cbetaselector,
    backbone_atom_names,
    backboneselector,
    heavyatomselector,
    resnameselector,
    water_res_names,
    waterselector,
    stdresselector,
    hetresselector,
    disorderselector,
    hydrogenselector,
    AminoAcidSequence


using Bio.Seq
import Bio.Seq.AminoAcidSequence


"A macromolecular structural element."
abstract StructuralElement


"""
An atom that is part of a macromolecule - either an `Atom` or a
`DisorderedAtom`.
"""
abstract AbstractAtom <: StructuralElement

"An atom that is part of a macromolecule."
immutable Atom <: AbstractAtom
    serial::Int
    name::String
    alt_loc_id::Char
    coords::Vector{Float64}
    occupancy::Float64
    temp_fac::Float64
    element::String
    charge::String
    residue::StructuralElement
end

"A container to hold different versions of the same atom."
immutable DisorderedAtom <: AbstractAtom
    alt_loc_ids::Dict{Char, Atom}
    default::Char
end


"""
A residue (amino acid) or other molecule - either a`Residue` or a
`DisorderedResidue`.
"""
abstract AbstractResidue <: StructuralElement

"A residue (amino acid) or other molecule."
immutable Residue <: AbstractResidue
    name::String
    number::Int
    ins_code::Char
    het_res::Bool # Does the residue consist of hetatoms?
    atom_list::Vector{String}
    atoms::Dict{String, AbstractAtom}
    chain::StructuralElement
end

"A container to hold different versions of the same residue (point mutations)."
immutable DisorderedResidue <: AbstractResidue
    names::Dict{String, Residue}
    default::String
end


"A chain (molecule) from a macromolecular structure."
immutable Chain <: StructuralElement
    id::Char
    res_list::Vector{String}
    residues::Dict{String, AbstractResidue}
    model::StructuralElement
end


"A conformation of a macromolecular structure."
immutable Model <: StructuralElement
    number::Int
    chains::Dict{Char, Chain}
    structure::StructuralElement
end


"A container for multiple `Model`s."
immutable ProteinStructure <: StructuralElement
    name::String
    models::Dict{Int, Model}
end


"""
A record for a single atom, e.g. as represented in a Protein Data Bank (PDB)
file
"""
immutable AtomRecord
    het_atom::Bool
    serial::Int
    atom_name::Compat.ASCIIString
    alt_loc_id::Char
    res_name::Compat.ASCIIString
    chain_id::Char
    res_number::Int
    ins_code::Char
    coords::Vector{Float64}
    occupancy::Float64
    temp_fac::Float64
    element::Compat.ASCIIString
    charge::Compat.ASCIIString
end


# Constructors without sub-elements

Residue(name::AbstractString, number::Integer, ins_code::Char, het_res::Bool, chain::Chain) = Residue(
    name, number, ins_code, het_res, [], Dict(), chain)

Chain(id::Char, model::Model) = Chain(id, [], Dict(), model)

Model(number::Integer, structure::ProteinStructure) = Model(number, Dict(), structure)
Model(structure::ProteinStructure) = Model(1, structure)

ProteinStructure(name::AbstractString) = ProteinStructure(name, Dict())
ProteinStructure() = ProteinStructure("")


"A `StructuralElement` or `Vector` of `StructuralElement`s."
typealias StructuralElementOrList Union{
        StructuralElement,
        Vector{AbstractAtom},
        Vector{Atom},
        Vector{DisorderedAtom},
        Vector{AbstractResidue},
        Vector{Residue},
        Vector{DisorderedResidue},
        Vector{Chain},
        Vector{Model}
    }


# Allow accessing sub elements contained in an element like a dictionary
# e.g. allows you to do res[atom_name] rather than res.atoms[atom_name]
# setindex! should be used with caution as it can lead to inconsistencies
# e.g. adding an atom to a residue atom dict does not update the atom list
# Hence it is used internally but not documented for general use

# Accessing a DisorderedAtom with a character returns the Atom with that alt
#   loc ID
Base.getindex(disordered_atom::DisorderedAtom, alt_loc_id::Char) = disordered_atom.alt_loc_ids[alt_loc_id]
function Base.setindex!(disordered_atom::DisorderedAtom, atom::Atom, alt_loc_id::Char)
    disordered_atom.alt_loc_ids[alt_loc_id] = atom
    return disordered_atom
end

function findatombyname(res::Residue, atom_name::AbstractString)
    if length(atom_name) == 4
        return res.atoms[atom_name]
    elseif length(atom_name) == 3
        if " $atom_name" in keys(res.atoms)
            return res.atoms[" $atom_name"]
        elseif "$atom_name " in keys(res.atoms)
            return res.atoms["$atom_name "]
        end
    elseif length(atom_name) == 2
        if " $atom_name " in keys(res.atoms)
            return res.atoms[" $atom_name "]
        elseif "  $atom_name" in keys(res.atoms)
            return res.atoms["  $atom_name"]
        elseif "$atom_name  " in keys(res.atoms)
            return res.atoms["$atom_name  "]
        end
    elseif length(atom_name) == 1
        if " $atom_name  " in keys(res.atoms)
            return res.atoms[" $atom_name  "]
        elseif "  $atom_name " in keys(res.atoms)
            return res.atoms["  $atom_name "]
        elseif "   $atom_name" in keys(res.atoms)
            return res.atoms["   $atom_name"]
        elseif "$atom_name   " in keys(res.atoms)
            return res.atoms["$atom_name   "]
        end
    end
    throw(KeyError(atom_name))
end

# Accessing a Residue with an AbstractString returns the AbstractAtom with that
#   atom name
Base.getindex(res::Residue, atom_name::AbstractString) = findatombyname(res, atom_name)
function Base.setindex!(res::Residue, atom::AbstractAtom, atom_name::AbstractString)
    res.atoms[atom_name] = atom
    return res
end

# Accessing a DisorderedResidue with an AbstractString returns the AbstractAtom
#   in the default Residue with that atom name
# This is not necessarily intuitive, it may be expected to return the Residue
#   with that residue name
# However this way accessing an AbstractResidue always returns an AbstractAtom
Base.getindex(disordered_res::DisorderedResidue, atom_name::AbstractString) = findatombyname(defaultresidue(disordered_res), atom_name)
function Base.setindex!(disordered_res::DisorderedResidue, atom::AbstractAtom, atom_name::AbstractString)
    defaultresidue(disordered_res)[atom_name] = atom
    return disordered_res
end

# Accessing a Chain with an AbstractString returns the AbstractResidue with that
#  residue ID
Base.getindex(chain::Chain, res_id::AbstractString) = chain.residues[res_id]
function Base.setindex!(chain::Chain, res::AbstractResidue, res_id::AbstractString)
    chain.residues[res_id] = res
    return chain
end

# Accessing a Chain with an Integer returns the AbstractResidue with that residue ID
#   converted to a string
Base.getindex(chain::Chain, res_number::Integer) = chain.residues[string(res_number)]
function Base.setindex!(chain::Chain, res::AbstractResidue, res_number::Integer)
    chain.residues[string(res_number)] = res
    return chain
end

# Accessing a Model with a Char returns the Chain with that chain ID
Base.getindex(model::Model, chain_id::Char) = model.chains[chain_id]
function Base.setindex!(model::Model, chain::Chain, chain_id::Char)
    model.chains[chain_id] = chain
    return model
end

# Accessing a ProteinStructure with an Integer returns the Model with that model
#   number
Base.getindex(struc::ProteinStructure, model_number::Integer) = struc.models[model_number]
function Base.setindex!(struc::ProteinStructure, model::Model, model_number::Integer)
    struc.models[model_number] = model
    return struc
end

# Accessing a ProteinStructure with a Char returns the Chain with that chain ID
#   on the default model
Base.getindex(struc::ProteinStructure, chain_id::Char) = defaultmodel(struc)[chain_id]
function Base.setindex!(struc::ProteinStructure, chain::Chain, chain_id::Char)
    defaultmodel(struc)[chain_id] = chain
    return struc
end


# Getters and setters for structural elements

"Get the serial number of an `AbstractAtom`."
serial(atom::Atom) = atom.serial
serial(disordered_atom::DisorderedAtom) = serial(defaultatom(disordered_atom))

"Get the atom name of an `AbstractAtom`, stripped of whitespace."
function atomname(atom::Atom; spaces::Bool=false)
    if spaces
        return atom.name
    else
        return strip(atom.name)
    end
end

atomname(disordered_atom::DisorderedAtom; spaces::Bool=false) = atomname(defaultatom(disordered_atom); spaces=spaces)

"Get the alternative location ID of an `AbstractAtom`."
altlocid(atom::Atom) = atom.alt_loc_id
altlocid(disordered_atom::DisorderedAtom) = defaultaltlocid(disordered_atom)

"Get the x coordinate of an `AbstractAtom`."
x(atom::Atom) = atom.coords[1]
x(disordered_atom::DisorderedAtom) = x(defaultatom(disordered_atom))

"""
Set the x coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
x!(atom::Atom, x::Real) = (atom.coords[1] = x; atom)
x!(disordered_atom::DisorderedAtom, x::Real) = x!(defaultatom(disordered_atom), x)

"Get the y coordinate of an `AbstractAtom`."
y(atom::Atom) = atom.coords[2]
y(disordered_atom::DisorderedAtom) = y(defaultatom(disordered_atom))

"""
Set the y coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
y!(atom::Atom, y::Real) = (atom.coords[2] = y; atom)
y!(disordered_atom::DisorderedAtom, y::Real) = y!(defaultatom(disordered_atom), y)

"Get the z coordinate of an `AbstractAtom`."
z(atom::Atom) = atom.coords[3]
z(disordered_atom::DisorderedAtom) = z(defaultatom(disordered_atom))

"""
Set the z coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
z!(atom::Atom, z::Real) = (atom.coords[3] = z; atom)
z!(disordered_atom::DisorderedAtom, z::Real) = z!(defaultatom(disordered_atom), z)

"Get the atomic coordinates of an `AbstractAtom` as a `Vector{Float64}`."
coords(atom::Atom) = atom.coords
coords(disordered_atom::DisorderedAtom) = coords(defaultatom(disordered_atom))

"""
Set the coordinates of an `AbstractAtom` to a `Vector{Float64}` of length 3. For
`DisorderedAtom`s only the default atom is updated.
"""
function coords!(atom::Atom, coords::Vector{Float64})
    @assert length(coords) == 3 "3 coordinates must be given"
    x!(atom, coords[1])
    y!(atom, coords[2])
    z!(atom, coords[3])
    return atom
end

coords!(disordered_atom::DisorderedAtom, coords::Vector{Float64}) = coords!(
    defaultatom(disordered_atom), coords)

"Get the occupancy of an `AbstractAtom`. Defaults to `1.0`."
occupancy(atom::Atom) = atom.occupancy
occupancy(disordered_atom::DisorderedAtom) = occupancy(defaultatom(disordered_atom))

"Get the temperature factor on an `AbstractAtom`. Defaults to `0.0`."
tempfac(atom::Atom) = atom.temp_fac
tempfac(disordered_atom::DisorderedAtom) = tempfac(defaultatom(disordered_atom))

"Get the element of an `AbstractAtom`. Defaults to `"  "`."
element(atom::Atom) = atom.element
element(disordered_atom::DisorderedAtom) = element(defaultatom(disordered_atom))

"Get the charge on an `AbstractAtom`. Defaults to `"  "`."
charge(atom::Atom) = atom.charge
charge(disordered_atom::DisorderedAtom) = charge(defaultatom(disordered_atom))

"Get the `Residue` that an `AbstractAtom` belongs to."
residue(atom::Atom) = atom.residue
residue(disordered_atom::DisorderedAtom) = residue(defaultatom(disordered_atom))

"""
Determines if an `AbstractAtom` represents a non-hetero atom, i.e. came from an
ATOM record.
"""
ishetatom(atom::Atom) = ishetres(residue(atom))
ishetatom(disordered_atom::DisorderedAtom) = ishetatom(defaultatom(disordered_atom))

"""
Determines if an `AbstractAtom` is a `DisorderedAtom`, i.e. if there are
multiple locations present for an atom.
"""
isdisorderedatom(::Atom) = false
isdisorderedatom(::DisorderedAtom) = true

"""
Get the alternative location ID of the default `Atom` in a `DisorderedAtom`. The
default is the highest occupancy, or lowest character alt loc ID for ties.
"""
defaultaltlocid(disordered_atom::DisorderedAtom) = disordered_atom.default

# Constructor acts as a setter for the default alt loc ID
function DisorderedAtom(disordered_atom::DisorderedAtom, default::Char)
    @assert default in altlocids(disordered_atom) "The new default alternative location ID must be present in the atom"
    return DisorderedAtom(disordered_atom.alt_loc_ids, default)
end

"""
Access the default `Atom` in a `DisorderedAtom`. The default is the highest
occupancy, or lowest alt loc ID character for ties.
"""
defaultatom(disordered_atom::DisorderedAtom) = disordered_atom[defaultaltlocid(disordered_atom)]

"""
Get the list of alternative location IDs in a `DisorderedAtom`, sorted by atom
serial.
"""
altlocids(disordered_atom::DisorderedAtom) = sort(collect(keys(
    disordered_atom.alt_loc_ids)), by= alt_loc_id -> disordered_atom[alt_loc_id])

"""
Get a descriptive atom ID for an `AbstractAtom` as a `Tuple` of the form
(full residue ID, residue name, atom name).
"""
atomid(atom::Atom) = (resid(atom; full=true), resname(atom), atomname(atom))
atomid(disordered_atom::DisorderedAtom) = atomid(defaultatom(disordered_atom))


"Get the residue name of an `AbstractAtom` or `AbstractResidue`."
resname(atom::Atom) = resname(residue(atom))
resname(disordered_atom::DisorderedAtom) = resname(defaultatom(disordered_atom))
resname(res::Residue) = res.name
resname(disordered_res::DisorderedResidue) = defaultresname(disordered_res)

"Get the residue number of an `AbstractAtom` or `AbstractResidue`."
resnumber(atom::Atom) = resnumber(residue(atom))
resnumber(disordered_atom::DisorderedAtom) = resnumber(defaultatom(disordered_atom))
resnumber(res::Residue) = res.number
resnumber(disordered_res::DisorderedResidue) = resnumber(defaultresidue(disordered_res))

"Get the insertion code of an `AbstractAtom` or `AbstractResidue`."
inscode(atom::Atom) = inscode(residue(atom))
inscode(disordered_atom::DisorderedAtom) = inscode(defaultatom(disordered_atom))
inscode(res::Residue) = res.ins_code
inscode(disordered_res::DisorderedResidue) = inscode(defaultresidue(disordered_res))

"""
Determines if an `AbstractResidue` represents a hetero molecule, i.e. consists
of HETATM records.
"""
ishetres(res::Residue) = res.het_res
ishetres(disordered_res::DisorderedResidue) = ishetres(defaultresidue(disordered_res))

"""
Determines if an `AbstractAtom` represents a hetero atom, i.e. came from a
HETATM record, or if an `AbstractResidue` represents a hetero molecule, i.e.
consists of HETATM records.
"""
ishetero(atom::AbstractAtom) = ishetatom(atom)
ishetero(res::AbstractResidue) = ishetres(res)

"Get the `Chain` that an `AbstractAtom` or `AbstractResidue` belongs to."
chain(atom::Atom) = chain(residue(atom))
chain(disordered_atom::DisorderedAtom) = chain(defaultatom(disordered_atom))
chain(res::Residue) = res.chain
chain(disordered_res::DisorderedResidue) = chain(defaultresidue(disordered_res))

"Get the chain ID of an `AbstractAtom`, `AbstractResidue` or `Chain`."
chainid(element::Union{AbstractResidue, AbstractAtom}) = chainid(chain(element))
chainid(chain::Chain) = chain.id















model(atom::Atom) = model(chain(atom))
model(disordered_atom::DisorderedAtom) = model(defaultatom(disordered_atom))
model(res::Residue) = model(chain(res))
model(disordered_res::DisorderedResidue) = model(defaultresidue(disordered_res))
model(chain::Chain) = chain.model

"""
Get the model number of a `Model`, `Chain`, `AbstractResidue` or `AbstractAtom`.
"""
modelnumber(model::Model) = model.number
modelnumber(element::Union{Chain, AbstractResidue, AbstractAtom}) = modelnumber(model(element))




structure(atom::Atom) = structure(model(atom))
structure(disordered_atom::DisorderedAtom) = structure(defaultatom(disordered_atom))
structure(res::Residue) = structure(model(res))
structure(disordered_res::DisorderedResidue) = structure(defaultresidue(disordered_res))
structure(chain::Chain) = structure(model(chain))
structure(model::Model) = model.structure

"""
Get the structure name of a `ProteinStructure`, or the `ProteinStructure`
associated with a `StructuralElement`.
"""
structurename(element::Union{Model, Chain, AbstractResidue, AbstractAtom}) = structurename(structure(element))
structurename(struc::ProteinStructure) = struc.name




"""
Get a descriptive residue ID for an `AbstractAtom` or `AbstractResidue`. Format
is residue number then insertion code with \"H_\" in front for hetero residues.
If `full` equals true the chain ID is also added after a colon. Examples are
\"50A\", \"H_20\" and \"10:A\".
"""
function resid(res::AbstractResidue; full::Bool=false)
    if ishetero(res)
        if full
            if inscode(res) == ' '
                return "H_$(resnumber(res)):$(chainid(res))"
            else
                return "H_$(resnumber(res))$(inscode(res)):$(chainid(res))"
            end
        else
            if inscode(res) == ' '
                return "H_$(resnumber(res))"
            else
                return "H_$(resnumber(res))$(inscode(res))"
            end
        end
    else
        if full
            if inscode(res) == ' '
                return "$(resnumber(res)):$(chainid(res))"
            else
                return "$(resnumber(res))$(inscode(res)):$(chainid(res))"
            end
        else
            if inscode(res) == ' '
                return "$(resnumber(res))"
            else
                return "$(resnumber(res))$(inscode(res))"
            end
        end
    end
end

resid(atom::AbstractAtom; full::Bool=false) = resid(residue(atom); full=full)

function resid(hetatm::Bool, resnum::Int, inscode::Char)
    if hetatm
        if inscode == ' '
            return "H_$resnum"
        else
            return "H_$resnum$inscode"
        end
    else
        if inscode == ' '
            return "$resnum"
        else
            return "$resnum$inscode"
        end
    end
end




"Get the sorted list of `AbstractAtom`s in an `AbstractResidue`."
atomnames(res::Residue) = res.atom_list

"Access the dictionary of `AbstractAtom`s in an `AbstractResidue`."
atoms(res::Residue) = res.atoms


"""
Determine if an `AbstractResidue` is a `DisorderedResidue`, i.e. there are
multiple residue names with the same residue ID.
"""
isdisorderedres(::Residue) = false
isdisorderedres(::DisorderedResidue) = true

"Access the `Residue` in a `DisorderedResidue` with a certain residue name."
disorderedres(disordered_res::DisorderedResidue, res_name::AbstractString) = disordered_res.names[res_name]

"""
Get the name of the default `Residue` in a `DisorderedResidue`. The default is
the first name read in.
"""
defaultresname(disordered_res::DisorderedResidue) = disordered_res.default

"""
Access the default `Residue` in a `DisorderedResidue`. The default is the first
name read in.
"""
defaultresidue(disordered_res::DisorderedResidue) = disordered_res.names[defaultresname(disordered_res)]

"""
Get the residue names in a `DisorderedResidue`. Default residue name is placed
first, then others are ordered alphabetically.
"""
resnames(disordered_res::DisorderedResidue) = sort(collect(keys(disordered_res.names)), lt= (res_name_one, res_name_two) -> (isless(res_name_one, res_name_two) && res_name_two != defaultresname(disordered_res)) || res_name_one == defaultresname(disordered_res))


atomnames(disordered_res::DisorderedResidue) = atomnames(defaultresidue(disordered_res))
atoms(disordered_res::DisorderedResidue) = atoms(defaultresidue(disordered_res))



# Constructor acts as a setter for the default residue name
function DisorderedResidue(disordered_res::DisorderedResidue, default::AbstractString)
    @assert default in resnames(disordered_res) "The new default residue name must be present in the residue"
    return DisorderedResidue(disordered_res.names, default)
end



"Get the sorted list of `AbstractResidue`s in a `Chain`."
resids(chain::Chain) = chain.res_list

"Access the dictionary of `AbstractResidue`s in a `Chain`."
residues(chain::Chain) = chain.residues




"Get the sorted chain IDs of the chains in a `Model` or `ProteinStructure`."
chainids(model::Model) = map(chainid, sort(collect(values(chains(model)))))

function chainids(struc::ProteinStructure)
    if countmodels(struc) > 0
        return chainids(defaultmodel(struc))
    else
        return Char[]
    end
end

"Access the dictionary of `Chain`s in a `Model`."
chains(model::Model) = model.chains


"Get the sorted `Model` numbers from a `ProteinStructure`."
modelnumbers(struc::ProteinStructure) = map(modelnumber, sort(collect(values(struc.models))))

"Access the dictionary of `Model`s in a `ProteinStructure`."
models(struc::ProteinStructure) = struc.models

"""
Get the default `Model` in a `ProteinStructure`. This is the `Model` with the
lowest model number.
"""
defaultmodel(struc::ProteinStructure) = struc.models[modelnumbers(struc)[1]]







# Sort lists of elements

# Sort atoms by serial
Base.isless(atom_one::AbstractAtom, atom_two::AbstractAtom) = isless(
    serial(atom_one), serial(atom_two))

# Sort residues by chain, then hetero, then resumber, then ins code
function Base.isless(res_one::AbstractResidue, res_two::AbstractResidue)
    if isless(chain(res_one), chain(res_two))
        return true
    elseif chainid(res_one) == chainid(res_two)
        if !ishetres(res_one) && ishetres(res_two)
            return true
        elseif ishetres(res_one) == ishetres(res_two)
            if isless(resnumber(res_one), resnumber(res_two))
                return true
            elseif resnumber(res_one) == resnumber(res_two)
                if isless(inscode(res_one), inscode(res_two))
                    return true
                end
            end
        end
    end
    return false
end

# Sort chains by chain ID
# Ordering is character sorting except the empty chain ID comes last
function Base.isless(chain_one::Chain, chain_two::Chain)
    if Int(chainid(chain_one)) < Int(chainid(chain_two)) && chainid(chain_one) != ' ' && chainid(chain_two) != ' '
        return true
    elseif chainid(chain_two) == ' ' && chainid(chain_one) != ' '
        return true
    else
        return false
    end
end

# Sort models by model number
Base.isless(model_one::Model, model_two::Model) = isless(
    modelnumber(model_one), modelnumber(model_two))


"""
Determine if the second residue follows the first in sequence. For this to be
true the residues need to be part of the same chain, both need to be
standard/hetero residues and the residue number of the second needs to be one
greater than that of the first.
"""
function sequentialresidues(res_first::AbstractResidue, res_second::AbstractResidue)
    if resnumber(res_second) == resnumber(res_first) + 1 && chain(res_second) == chain(res_first) && ishetres(res_second) == ishetres(res_first)
        return true
    else
        return false
    end
end


# Iterators to yield sub elements when looping over an element

# Iterating over a ProteinStructure yields Models
Base.length(struc::ProteinStructure) = length(modelnumbers(struc))
Base.start(::ProteinStructure) = 1
Base.next(struc::ProteinStructure, state) = (struc[modelnumbers(struc)[state]], state + 1)
Base.done(struc::ProteinStructure, state) = state > length(struc)
Base.eltype(::Type{ProteinStructure}) = Model

# Iterating over a Model yields Chains
Base.length(model::Model) = length(chainids(model))
Base.start(::Model) = 1
Base.next(model::Model, state) = (model[chainids(model)[state]], state + 1)
Base.done(model::Model, state) = state > length(model)
Base.eltype(::Type{Model}) = Chain

# Iterating over a Chain yields AbstractResidues
Base.length(chain::Chain) = length(resids(chain))
Base.start(::Chain) = 1
Base.next(chain::Chain, state) = (chain[resids(chain)[state]], state + 1)
Base.done(chain::Chain, state) = state > length(chain)
Base.eltype(::Type{Chain}) = AbstractResidue

# Iterating over a Residue yields AbstractAtoms
Base.length(res::Residue) = length(atomnames(res))
Base.start(::Residue) = 1
Base.next(res::Residue, state) = (res[atomnames(res)[state]], state + 1)
Base.done(res::Residue, state) = state > length(res)
Base.eltype(::Type{Residue}) = AbstractAtom

# Iterating over a DisorderedResidue yields AbstractAtoms
# This is not necessarily intuitive, it may be expected to yield Residues
# However this way iterating over an AbstractResidue always yields AbstractAtoms
Base.length(disordered_res::DisorderedResidue) = length(atomnames(disordered_res))
Base.start(::DisorderedResidue) = 1
Base.next(disordered_res::DisorderedResidue, state) = (defaultresidue(disordered_res)[atomnames(disordered_res)[state]], state + 1)
Base.done(disordered_res::DisorderedResidue, state) = state > length(disordered_res)
Base.eltype(::Type{DisorderedResidue}) = AbstractAtom

# Iterating over an Atom returns itself
# This is not necessarily intuitive, it may be expected to not be an iterator
# However this way iterating over an AbstractAtom always yields Atoms
Base.length(::Atom) = 1
Base.start(::Atom) = 1
Base.next(atom::Atom, state) = (atom, state + 1)
Base.done(::Atom, state) = state > 1
Base.eltype(::Type{Atom}) = Atom

# Iterating over a DisorderedAtom yields Atoms
Base.length(disordered_atom::DisorderedAtom) = length(altlocids(disordered_atom))
Base.start(::DisorderedAtom) = 1
Base.next(disordered_atom::DisorderedAtom, state) = (disordered_atom[altlocids(disordered_atom)[state]], state + 1)
Base.done(disordered_atom::DisorderedAtom, state) = state > length(disordered_atom)
Base.eltype(::Type{DisorderedAtom}) = Atom


"""
Returns a copy of a `Vector` of atoms or residues with all elements that do not
satisfy `selector_functions...` removed.
"""
function applyselectors{T <: Union{AbstractResidue, AbstractAtom}}(element_list::Vector{T}, selector_functions::Function...)
    new_list = copy(element_list)
    applyselectors!(new_list, selector_functions...)
    return new_list
end

"Runs `applyselectors` in place."
function applyselectors!{T <: Union{AbstractResidue, AbstractAtom}}(element_list::Vector{T}, selector_functions::Function...)
    for selector_function in selector_functions
        filter!(selector_function, element_list)
    end
    return element_list
end


"""
Returns a sorted `Vector` of the residues in a `StructuralElementOrList`.
Additional arguments are `selector_functions...` - only residues that satisfy
the selector functions are retained.
"""
function collectresidues(struc::ProteinStructure, selector_functions::Function...)
    if countmodels(struc) > 0
        return collectresidues(defaultmodel(struc), selector_functions...)
    else
        return AbstractResidue[]
    end
end

collectresidues(chain::Chain, selector_functions::Function...) = applyselectors(collect(chain), selector_functions...)
collectresidues(res::AbstractResidue, selector_functions::Function...) = applyselectors(AbstractResidue[res], selector_functions...)
collectresidues(atom::AbstractAtom, selector_functions::Function...) = applyselectors(organise(atom), selector_functions...)
# Note output is always Vector{AbstractResidue} unless input was Vector{Residue}
# or Vector{DisorderedResidue}, in which case output is same type as input type
collectresidues{T <: AbstractResidue}(residues::Vector{T}, selector_functions::Function...) = sort(applyselectors(residues, selector_functions...))
collectresidues{T <: AbstractAtom}(atoms::Vector{T}, selector_functions::Function...) = applyselectors(organise(atoms), selector_functions...)

function collectresidues(element::Union{Model, Vector{Model}, Vector{Chain}}, selector_functions::Function...)
    residues = AbstractResidue[]
    for sub_element in element
        append!(residues, collectresidues(sub_element, selector_functions...))
    end
    return residues
end


"""
Returns a sorted `Vector` of the atoms in a `StructuralElementOrList`.
Additional arguments are `selector_functions...` - only atoms that satisfy the
selector functions are retained.
"""
function collectatoms(struc::ProteinStructure, selector_functions::Function...)
    if countmodels(struc) > 0
        return collectatoms(defaultmodel(struc), selector_functions...)
    else
        return AbstractAtom[]
    end
end

collectatoms(res::AbstractResidue, selector_functions::Function...) = applyselectors(collect(res), selector_functions...)
collectatoms(atom::AbstractAtom, selector_functions::Function...) = applyselectors(AbstractAtom[atom], selector_functions...)
# Note output is always Vector{AbstractAtom} unless input was Vector{Atom} or
# Vector{DisorderedAtom}, in which case output is same type as input type
collectatoms{T <: AbstractAtom}(atoms::Vector{T}, selector_functions::Function...) = sort(applyselectors(atoms, selector_functions...))

function collectatoms(element::Union{Model, Chain, Vector{Model}, Vector{Chain}, Vector{AbstractResidue}, Vector{Residue}, Vector{DisorderedResidue}}, selector_functions::Function...)
    atoms = AbstractAtom[]
    for sub_element in element
        append!(atoms, collectatoms(sub_element, selector_functions...))
    end
    return atoms
end


"Get the number of `Model`s in a `ProteinStructure`."
countmodels(struc::ProteinStructure) = length(struc)

"Get the number of `Chain`s in a `ProteinStructure` or a `Model`."
function countchains(struc::ProteinStructure)
    if countmodels(struc) > 0
        return length(defaultmodel(struc))
    else
        return 0
    end
end

countchains(model::Model) = length(model)

"""
Get the number of residues in a `StructuralElementOrList`. Additional arguments
are `selector_functions...` - only residues that satisfy the selector functions
are counted.
"""
countresidues(element::StructuralElementOrList, selector_functions::Function...) = length(collectresidues(element, selector_functions...))

"""
Get the number of atoms in a `StructuralElementOrList`. Additional arguments are `selector_functions...` - only atoms that satisfy the selector functions are
counted.
"""
countatoms(element::StructuralElementOrList, selector_functions::Function...) = length(collectatoms(element, selector_functions...))


"""
Organise a `StructuralElementOrList` into the next level up the heirarchy. A
`Vector{AbstractAtom}` becomes a `Vector{AbstractResidue}`, a
`Vector{AbstractResidue}` becomes a `Vector{Chain}`, a `Vector{Chain}` becomes a
`Vector{Model}` and a `Vector{Model}` becomes a `Vector{ProteinStructure}`.
"""
function organise(models::Vector{Model})
    strucs_out = ProteinStructure[]
    for model in models
        if !(structure(model) in strucs_out)
            push!(strucs_out, structure(model))
        end
    end
    return strucs_out
end

function organise(chains::Vector{Chain})
    models_out = Model[]
    for chain in chains
        if !(model(chain) in models_out)
            push!(models_out, model(chain))
        end
    end
    return sort(models_out)
end

function organise{T <: AbstractResidue}(residues::Vector{T})
    chains_out = Chain[]
    for res in residues
        if !(chain(res) in chains_out)
            push!(chains_out, chain(res))
        end
    end
    return sort(chains_out)
end

function organise{T <: AbstractAtom}(atoms::Vector{T})
    residues_out = AbstractResidue[]
    for atom in atoms
        if !(residue(atom) in residues_out)
            push!(residues_out, residue(atom))
        end
    end
    return sort(residues_out)
end

organise(element::Union{Model, Chain, AbstractResidue, AbstractAtom}) = organise([element])

# Throw errors when struc passed?
function organisemodel{T <: StructuralElement}(element_list::Vector{T})
    models_out = Model[]
    for element in element_list
        if !(model(element) in models_out)
            push!(models_out, model(element))
        end
    end
    return sort(models_out)
end

organisemodel(element::StructuralElement) = organisemodel([element])

function organisestructure{T <: StructuralElement}(element_list::Vector{T})
    strucs_out = ProteinStructure[]
    for element in element_list
        if !(structure(element) in strucs_out)
            push!(strucs_out, structure(element))
        end
    end
    return strucs_out
end

organisestructure(element::StructuralElement) = organisestructure([element])


function unsafe_addatomtomodel!(model::Model, atom_record::AtomRecord; remove_disorder::Bool=false)
    # Add chain to model if necessary
    if !haskey(chains(model), atom_record.chain_id)
        model[atom_record.chain_id] = Chain(atom_record.chain_id, model)
    end
    chain = model[atom_record.chain_id]
    res_id = resid(atom_record.het_atom, atom_record.res_number, atom_record.ins_code)
    # If residue does not exist in the chain, create a Residue
    if !haskey(residues(chain), res_id)
        chain[res_id] = Residue(atom_record.res_name, atom_record.res_number, atom_record.ins_code, atom_record.het_atom, chain)
        residue = chain[res_id]
    elseif isa(chain[res_id], Residue)
        # Residue exists in the chain and the residue names match
        # Add to that Residue
        if resname(chain[res_id]) == atom_record.res_name
            residue = chain[res_id]
        # Residue exists in the chain but the residue names do not match
        # Create a DisorderedResidue
        else
            chain[res_id] = DisorderedResidue(Dict(
                resname(chain[res_id]) => chain[res_id],
                atom_record.res_name => Residue(atom_record.res_name, atom_record.res_number, atom_record.ins_code, atom_record.het_atom, chain)
            ), resname(chain[res_id]))
            residue = disorderedres(chain[res_id], atom_record.res_name)
        end
    else
        # DisorderedResidue exists in the chain and the residue names match
        # Add to that DisorderedResidue
        if atom_record.res_name in resnames(chain[res_id])
            residue = disorderedres(chain[res_id], atom_record.res_name)
        # DisorderedResidue exists in the chain and the residue names do not match
        # Create a new Residue in the DisorderedResidue
        else
            chain[res_id].names[atom_record.res_name] = Residue(atom_record.res_name, atom_record.res_number, atom_record.ins_code, atom_record.het_atom, chain) # Make residues function for this?
            residue = disorderedres(chain[res_id], atom_record.res_name)
        end
    end
    atom = Atom(atom_record.serial, atom_record.atom_name, atom_record.alt_loc_id, atom_record.coords, atom_record.occupancy, atom_record.temp_fac, atom_record.element, atom_record.charge, residue)
    # If atom does not exist in the residue, create an Atom
    if !haskey(atoms(residue), atom_record.atom_name)
        residue[atom_record.atom_name] = atom
    # Atom exists in the residue, atom names match and alt loc IDs are different
    elseif isa(residue[atom_record.atom_name], Atom) && atom_record.alt_loc_id != altlocid(residue[atom_record.atom_name])
        # If we are removing disorder and the new atom is preferred to the old one, replace the old one
        if remove_disorder && choosedefaultaltlocid(atom, residue[atom_record.atom_name]) == atom_record.alt_loc_id
             residue[atom_record.atom_name] = atom
        # If we are not removing disorder, create a new disordered atom container and add both atoms
        elseif !remove_disorder
            residue[atom_record.atom_name] = DisorderedAtom(Dict(
                atom_record.alt_loc_id => atom,
                altlocid(residue[atom_record.atom_name]) => residue[atom_record.atom_name]
            ), choosedefaultaltlocid(atom, residue[atom_record.atom_name]))
        end
    # A disordered atom container already exists and the alt loc ID is not taken
    elseif isa(residue[atom_record.atom_name], DisorderedAtom) && !(atom_record.alt_loc_id in altlocids(residue[atom_record.atom_name]))
        # Add the new atom to the disordered atom container
        residue[atom_record.atom_name][atom_record.alt_loc_id] = atom
        # If the default alt loc requires changing, change it
        if choosedefaultaltlocid(defaultatom(residue[atom_record.atom_name]), atom) != defaultaltlocid(residue[atom_record.atom_name])
            residue[atom_record.atom_name] = DisorderedAtom(residue[atom_record.atom_name], atom_record.alt_loc_id)
        end
    else
        error("Two copies of the same atom have the same alternative location ID:\n$(residue[atom_record.atom_name])\n$(atom)")
    end
end


# temp name
function fixlists!(struc::ProteinStructure)
    for model in struc
        for chain in model
            append!(chain.res_list, map(resid, sort(collect(values(residues(chain))))))
            for res in chain
                if isa(res, Residue)
                    fixlists!(res)
                else
                    for res_name in resnames(res)
                        fixlists!(disorderedres(res, res_name))
                    end
                end
            end
        end
    end
end


fixlists!(res::Residue) = append!(res.atom_list, map(atomname, sort(collect(values(atoms(res))))))


"""
Determine which of two `Atom`s representing a disorered atom better qualifies
as the default location. The `Atom` with the highest occupancy is chosen; in the
case of ties the `Atom` with the lowest alternative location ID in alphabetical
order is chosen.
"""
function choosedefaultaltlocid(atom_one::Atom, atom_two::Atom)
    if occupancy(atom_one) > occupancy(atom_two) ||
            (occupancy(atom_one) == occupancy(atom_two) &&
            Int(altlocid(atom_one)) < Int(altlocid(atom_two)))
        return altlocid(atom_one)
    else
        return altlocid(atom_two)
    end
end


"""
Determines if an `AbstractAtom` represents a non-hetero atom, i.e. came from an
ATOM record.
"""
stdatomselector(atom::AbstractAtom) = !ishetatom(atom)

"""
Determines if an `AbstractAtom` is a hetero atom, i.e. came from a HETATM
record.
"""
hetatomselector(atom::AbstractAtom) = ishetatom(atom)

"""
Determines if an `AbstractAtom` has its atom name in the given `Set` or
`Vector`.
"""
atomnameselector(atom::AbstractAtom, atom_names::Set{String}) = atomname(atom) in atom_names
# Set is faster but Vector method is retained for ease of use
atomnameselector(atom::AbstractAtom, atom_names::Vector{String}) = atomname(atom) in atom_names

"`Set` of C-alpha atom names."
const calpha_atom_names = Set(["CA"])

"""
Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
C-alpha atom.
"""
calphaselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, calpha_atom_names)

"`Set` of C-beta atom names."
const cbeta_atom_names = Set(["CB"])

"""
Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
C-beta atom, or a C-alpha atom in glycine.
"""
cbetaselector(atom::AbstractAtom) = stdatomselector(atom) && (atomnameselector(atom, cbeta_atom_names) || (atomnameselector(atom, calpha_atom_names) && resname(atom) == "GLY"))

"`Set` of protein backbone atom names."
const backbone_atom_names = Set(["CA", "N", "C"])

"""
Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
protein backbone atom.
"""
backboneselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, backbone_atom_names)

"""
Determines if an `AbstractAtom` corresponds to a heavy (non-hydrogen) atom and
is not a hetero-atom.
"""
heavyatomselector(atom::AbstractAtom) = stdatomselector(atom) && !hydrogenselector(atom)

"""
Determines if an `AbstractResidue` or `AbstractAtom` has its resiudue name in
the given `Set` or `Vector`.
"""
resnameselector(element::Union{AbstractResidue, AbstractAtom}, res_names::Set{String}) = resname(element) in res_names
resnameselector(element::Union{AbstractResidue, AbstractAtom}, res_names::Vector{String}) = resname(element) in res_names

"`Set` of residue names corresponding to water."
const water_res_names = Set(["HOH"])

"""
Determines if an `AbstractResidue` or `AbstractAtom` represents a water
molecule.
"""
waterselector(element::Union{AbstractResidue, AbstractAtom}) = resnameselector(element, water_res_names)

"""
Determines if an `AbstractResidue` represents a standard protein residue,
i.e. consists of ATOM records.
"""
stdresselector(res::AbstractResidue) = !ishetres(res)

"""
Determines if an `AbstractResidue` represents a hetero molecule, i.e. consists
of HETATM records.
"""
hetresselector(res::AbstractResidue) = ishetres(res)

"""
Determines whether an `AbstractAtom` or `AbstractResidue` is disordered, i.e.
has multiple locations in the case of atoms or multiple residue names (point
mutants) in the case of residues.
"""
disorderselector(atom::AbstractAtom) = isdisorderedatom(atom)
disorderselector(res::AbstractResidue) = isdisorderedres(res)

# Either the element is H or the element field is empty, the atom name contains
#   an H and there are no letters in the atom name before H
# For example atom names 1H and H1 would be hydrogens but NH1 would not
"""
Determines if an `AbstractAtom` represents hydrogen. Uses the element field
where possible, otherwise uses the atom name.
"""
hydrogenselector(atom::AbstractAtom) = element(atom) == "H" || (element(atom) == "" && 'H' in atomname(atom) && !ismatch(r"[a-zA-Z]", atomname(atom)[1:findfirst(atomname(atom), 'H')-1]))


# how to deal with a missing res? ignore or gap?
# what should unknown res be?
# have file as input and wrap?
function AminoAcidSequence(element::StructuralElementOrList, selector_functions::Function...)
    residues = collectresidues(element, selector_functions...)
    seq = AminoAcid[]
    for res in residues
        try
            push!(seq, parse(AminoAcid, resname(res)))
        catch
            push!(seq, AA_X)
        end
    end
    return AminoAcidSequence(seq)
end


# Descriptive showing of elements

function Base.show(io::IO, struc::ProteinStructure)
    if countmodels(struc) > 0
        model = defaultmodel(struc)
        println(io, "Name                        -  ", structurename(struc))
        println(io, "Number of models            -  ", countmodels(struc))
        println(io, "Chain(s)                    -  ", join(chainids(model)))
        println(io, "Number of residues          -  ", countresidues(model, stdresselector))
        println(io, "Number of point mutations   -  ", countresidues(model, stdresselector, disorderselector))
        println(io, "Number of other molecules   -  ", countresidues(model, hetresselector) - countresidues(model, hetresselector, waterselector))
        println(io, "Number of water molecules   -  ", countresidues(model, hetresselector, waterselector))
        println(io, "Number of atoms             -  ", countatoms(model, stdatomselector))
        println(io, "Number of hydrogens         -  ", countatoms(model, stdatomselector, hydrogenselector))
        println(io, "Number of disordered atoms  -  ", countatoms(model, stdatomselector, disorderselector))
    else
        println(io, "Name                        -  ", structurename(struc))
        println(io, "Number of models            -  0")
    end
end

function Base.show(io::IO, model::Model)
    println(io, "Model number                -  ", modelnumber(model))
    println(io, "Chain(s)                    -  ", join(chainids(model)))
    println(io, "Number of residues          -  ", countresidues(model, stdresselector))
    println(io, "Number of point mutations   -  ", countresidues(model, stdresselector, disorderselector))
    println(io, "Number of other molecules   -  ", countresidues(model, hetresselector) - countresidues(model, hetresselector, waterselector))
    println(io, "Number of water molecules   -  ", countresidues(model, hetresselector, waterselector))
    println(io, "Number of atoms             -  ", countatoms(model, stdatomselector))
    println(io, "Number of hydrogens         -  ", countatoms(model, stdatomselector, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(model, stdatomselector, disorderselector))
end

function Base.show(io::IO, chain::Chain)
    println(io, "Chain ID                    -  ", chainid(chain))
    println(io, "Number of residues          -  ", countresidues(chain, stdresselector))
    println(io, "Number of point mutations   -  ", countresidues(chain, stdresselector, disorderselector))
    println(io, "Number of other molecules   -  ", countresidues(chain, hetresselector) - countresidues(chain, hetresselector, waterselector))
    println(io, "Number of water molecules   -  ", countresidues(chain, hetresselector, waterselector))
    println(io, "Number of atoms             -  ", countatoms(chain, stdatomselector))
    println(io, "Number of hydrogens         -  ", countatoms(chain, stdatomselector, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(chain, stdatomselector, disorderselector))
end

function Base.show(io::IO, res::Residue)
    println(io, "Residue ID                  -  ", resid(res; full=true))
    println(io, "Residue name                -  ", resname(res))
    println(io, "Number of atoms             -  ", countatoms(res))
    println(io, "Number of hydrogens         -  ", countatoms(res, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(res, disorderselector))
end

function Base.show(io::IO, disordered_res::DisorderedResidue)
    println(io, "Residue ID                  -  ", resid(disordered_res; full=true))
    for res_name in resnames(disordered_res)
        println(io, "Residue name                -  ", res_name)
        println(io, "Number of atoms             -  ", countatoms(disorderedres(disordered_res, res_name)))
        println(io, "Number of hydrogens         -  ", countatoms(disorderedres(disordered_res, res_name), hydrogenselector))
        println(io, "Number of disordered atoms  -  ", countatoms(disorderedres(disordered_res, res_name), disorderselector))
    end
end

Base.show(io::IO, atom::Atom) = println(io, pdbline(atom)...)
Base.showcompact(io::IO, atom::Atom) = print(io, pdbline(atom)...)

function Base.show(io::IO, disordered_atom::DisorderedAtom)
    for atom in disordered_atom
        show(io, atom)
    end
end
