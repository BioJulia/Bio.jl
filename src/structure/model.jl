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
    sortchainids,
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


"A protein structural element."
abstract StructuralElement


"An atom represented in a PDB file - either an `Atom` or a `DisorderedAtom`."
abstract AbstractAtom <: StructuralElement

"An atom represented by an ATOM or HETATM record from a PDB file."
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


"A container to hold different versions of the same atom in a PDB file."
immutable DisorderedAtom <: AbstractAtom
    alt_loc_ids::Dict{Char, Atom}
    default::Char
end


"""
A residue (amino acid) or other molecule represented in a PDB file - either a
`Residue` or a `DisorderedResidue`.
"""
abstract AbstractResidue <: StructuralElement

"A residue (amino acid) or other molecule record from a PDB file."
immutable Residue <: AbstractResidue
    name::String
    number::Int
    ins_code::Char
    het_res::Bool # Does the residue consist of hetatoms?
    atom_list::Vector{String}
    atoms::Dict{String, AbstractAtom}
    chain::StructuralElement
end

# Constructor from a list of atoms
Residue{T <: AbstractAtom}(atoms::Vector{T}) = Residue(
        resname(atoms[1]),
        chainid(atoms[1]),
        resnumber(atoms[1]),
        inscode(atoms[1]),
        ishetero(atoms[1]),
        map(atomname, sort(atoms)),
        Dict(atomname(atom) => atom for atom in atoms))

"A container to hold different versions of the same residue (point mutations)."
immutable DisorderedResidue <: AbstractResidue
    names::Dict{String, Residue}
    default::String
end


"A chain represented in a PDB file."
immutable Chain <: StructuralElement
    id::Char
    res_list::Vector{String}
    residues::Dict{String, AbstractResidue}
    model::StructuralElement
end


"A conformation of a protein represented in a PDB file."
immutable Model <: StructuralElement
    number::Int
    chains::Dict{Char, Chain}
    structure::StructuralElement
end


"A container for multiple models from a PDB file."
immutable ProteinStructure <: StructuralElement
    name::String
    models::Dict{Int, Model}
end

ProteinStructure(name::AbstractString) = ProteinStructure(name, Dict())
ProteinStructure() = ProteinStructure("")

# Moved so types defined in the right order
# Constructor without atoms
Residue(name::AbstractString, number::Integer, ins_code::Char, het_res::Bool, chain::Chain) = Residue(
    name, number, ins_code, het_res, [], Dict(), chain)

Chain(id::Char, model::Model) = Chain(id, [], Dict(), model)

Model(number::Integer, structure::ProteinStructure) = Model(number, Dict(), structure)
Model(structure::ProteinStructure) = Model(1, structure)


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

# Accessing a DisorderedAtom with a character returns the Atom with that alt
#   loc ID
Base.getindex(disordered_atom::DisorderedAtom, alt_loc_id::Char) = disordered_atom.alt_loc_ids[alt_loc_id]
function Base.setindex!(disordered_atom::DisorderedAtom, atom::Atom, alt_loc_id::Char)
    disordered_atom.alt_loc_ids[alt_loc_id] = atom
    return disordered_atom
end

# Accessing a Residue with an AbstractString returns the AbstractAtom with that
#   atom name
Base.getindex(res::Residue, atom_name::AbstractString) = res.atoms[atom_name]
function Base.setindex!(res::Residue, atom::AbstractAtom, atom_name::AbstractString)
    res.atoms[atom_name] = atom
    return res
end

# Accessing a DisorderedResidue with an AbstractString returns the AbstractAtom
#   in the default Residue with that atom name
# This is not necessarily intuitive, it may be expected to return the Residue
#   with that residue name
# However this way accessing an AbstractResidue always returns an AbstractAtom
Base.getindex(disordered_res::DisorderedResidue, atom_name::AbstractString) = disordered_res.names[defaultresname(disordered_res)][atom_name]
function Base.setindex!(disordered_res::DisorderedResidue, atom::AbstractAtom, atom_name::AbstractString)
    disordered_res.names[defaultresname(disordered_res)][atom_name] = atom
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

# Accessing a ProteinStructure with an Integer returns the Model with that model
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

# Atom getters/setters

"""
Determines if an `AbstractAtom` represents a non-hetero atom, i.e. came from an
ATOM record.
"""
ishetatom(atom::Atom) = ishetres(residue(atom))

"Get the serial number of an `AbstractAtom`."
serial(atom::Atom) = atom.serial

"Get the atom name of an `AbstractAtom`, stripped of whitespace."
atomname(atom::Atom) = atom.name

"Get the alternative location ID of an `AbstractAtom`."
altlocid(atom::Atom) = atom.alt_loc_id

"Get the residue name of an `AbstractAtom` or `AbstractResidue`."
resname(atom::Atom) = resname(residue(atom))

"Get the chain ID of an `AbstractAtom`, `AbstractResidue` or `Chain`."
chainid(element::Union{AbstractResidue, AbstractAtom}) = chainid(chain(element))

"Get the residue number of an `AbstractAtom` or `AbstractResidue`."
resnumber(atom::Atom) = resnumber(residue(atom))

"Get the insertion code of an `AbstractAtom` or `AbstractResidue`."
inscode(atom::Atom) = inscode(residue(atom))

"Get the x coordinate of an `AbstractAtom`."
x(atom::Atom) = atom.coords[1]

"Get the y coordinate of an `AbstractAtom`."
y(atom::Atom) = atom.coords[2]

"Get the z coordinate of an `AbstractAtom`."
z(atom::Atom) = atom.coords[3]

"Get the atomic coordinates of an `AbstractAtom` as a `Vector{Float64}`."
coords(atom::Atom) = atom.coords

"""
Get the occupancy of an `AbstractAtom`. Defaults to `1.0` if not read from the
PDB file.
"""
occupancy(atom::Atom) = atom.occupancy

"""
Get the temperature factor on an `AbstractAtom`. Defaults to `0.0` if not read
from the PDB file.
"""
tempfac(atom::Atom) = atom.temp_fac

"""
Get the element of an `AbstractAtom`. Defaults to `""` if not read from the PDB
file.
"""
element(atom::Atom) = atom.element

"""
Get the charge on an `AbstractAtom`. Defaults to `""` if not read from the PDB
file.
"""
charge(atom::Atom) = atom.charge

residue(atom::Atom) = atom.residue
chain(atom::Atom) = chain(residue(atom))
model(atom::Atom) = model(chain(atom))
structure(atom::Atom) = structure(model(atom))

modelnumber(element::Union{Chain, AbstractResidue, AbstractAtom}) = modelnumber(model(element))
structurename(element::Union{Model, Chain, AbstractResidue, AbstractAtom}) = structurename(structure(element))

"""
Determines if an `AbstractAtom` represents a hetero atom, i.e. came from a
HETATM record, or if an `AbstractResidue` represents a hetero molecule, i.e.
consists of HETATM records.
"""
ishetero(atom::AbstractAtom) = ishetatom(atom)

"""
Determines if an `AbstractAtom` is a `DisorderedAtom`, i.e. if there are
multiple locations present for an atom.
"""
isdisorderedatom(::Atom) = false

"""
Get a descriptive atom ID for an `AbstractAtom` as a `Tuple` of the form
(full residue ID, residue name, atom name).
"""
atomid(atom::Atom) = (resid(atom; full=true), resname(atom), atomname(atom))


"""
Set the x coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
x!(atom::Atom, x::Real) = (atom.coords[1] = x; atom)

"""
Set the y coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
y!(atom::Atom, y::Real) = (atom.coords[2] = y; atom)

"""
Set the z coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
z!(atom::Atom, z::Real) = (atom.coords[3] = z; atom)


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

# DisorderedAtom getters/setters

"""
Get the alternative location ID of the default `Atom` in a `DisorderedAtom`. The
default is the highest occupancy, or lowest character alt loc ID for ties.
"""
defaultaltlocid(disordered_atom::DisorderedAtom) = disordered_atom.default

"""
Access the default `Atom` in a `DisorderedAtom`. The default is the highest
occupancy, or lowest alt loc ID character for ties.
"""
defaultatom(disordered_atom::DisorderedAtom) = disordered_atom[defaultaltlocid(disordered_atom)]

"""
Get the list of alternative location IDs in a `DisorderedAtom`, sorted by atom
serial.
"""
altlocids(disordered_atom::DisorderedAtom) = sort(collect(keys(disordered_atom.alt_loc_ids)), by= alt_loc_id -> serial(disordered_atom[alt_loc_id]))

ishetatom(disordered_atom::DisorderedAtom) = ishetatom(defaultatom(disordered_atom))
serial(disordered_atom::DisorderedAtom) = serial(defaultatom(disordered_atom))
atomname(disordered_atom::DisorderedAtom) = atomname(defaultatom(disordered_atom))
altlocid(disordered_atom::DisorderedAtom) = defaultaltlocid(disordered_atom)
resname(disordered_atom::DisorderedAtom) = resname(defaultatom(disordered_atom))
resnumber(disordered_atom::DisorderedAtom) = resnumber(defaultatom(disordered_atom))
inscode(disordered_atom::DisorderedAtom) = inscode(defaultatom(disordered_atom))
x(disordered_atom::DisorderedAtom) = x(defaultatom(disordered_atom))
y(disordered_atom::DisorderedAtom) = y(defaultatom(disordered_atom))
z(disordered_atom::DisorderedAtom) = z(defaultatom(disordered_atom))
coords(disordered_atom::DisorderedAtom) = coords(defaultatom(disordered_atom))
occupancy(disordered_atom::DisorderedAtom) = occupancy(defaultatom(disordered_atom))
tempfac(disordered_atom::DisorderedAtom) = tempfac(defaultatom(disordered_atom))
element(disordered_atom::DisorderedAtom) = element(defaultatom(disordered_atom))
charge(disordered_atom::DisorderedAtom) = charge(defaultatom(disordered_atom))

residue(disordered_atom::DisorderedAtom) = residue(defaultatom(disordered_atom))
chain(disordered_atom::DisorderedAtom) = chain(defaultatom(disordered_atom))
model(disordered_atom::DisorderedAtom) = model(defaultatom(disordered_atom))
structure(disordered_atom::DisorderedAtom) = structure(defaultatom(disordered_atom))

isdisorderedatom(::DisorderedAtom) = true
atomid(disordered_atom::DisorderedAtom) = atomid(defaultatom(disordered_atom))

# Constructor acts as a setter for the default alt loc ID
function DisorderedAtom(disordered_atom::DisorderedAtom, default::Char)
    @assert default in altlocids(disordered_atom) "The new default alternative location ID must be present in the atom"
    return DisorderedAtom(disordered_atom.alt_loc_ids, default)
end

# These coordinate setters only set the default atom coordinates
x!(disordered_atom::DisorderedAtom, x::Real) = x!(defaultatom(disordered_atom), x)
y!(disordered_atom::DisorderedAtom, y::Real) = y!(defaultatom(disordered_atom), y)
z!(disordered_atom::DisorderedAtom, z::Real) = z!(defaultatom(disordered_atom), z)
coords!(disordered_atom::DisorderedAtom, coords::Vector{Float64}) = coords!(defaultatom(disordered_atom), coords)

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


# Residue getters/setters
resname(res::Residue) = res.name
resnumber(res::Residue) = res.number
inscode(res::Residue) = res.ins_code

"""
Determines if an `AbstractResidue` represents a hetero molecule, i.e. consists
of HETATM records.
"""
ishetres(res::Residue) = res.het_res

chain(res::Residue) = res.chain
model(res::Residue) = model(chain(res))
structure(res::Residue) = structure(model(res))

"Get the sorted list of `AbstractAtom`s in an `AbstractResidue`."
atomnames(res::Residue) = res.atom_list

"Access the dictionary of `AbstractAtom`s in an `AbstractResidue`."
atoms(res::Residue) = res.atoms

ishetero(res::AbstractResidue) = ishetres(res)

"""
Determine if an `AbstractResidue` is a `DisorderedResidue`, i.e. there are
multiple residue names with the same residue ID.
"""
isdisorderedres(::Residue) = false

# DisorderedResidue getters/setters

"Access the `Residue` in a `DisorderedResidue` with a certain residue name."
disorderedres(disordered_res::DisorderedResidue, res_name::AbstractString) = disordered_res.names[res_name]

"""
Get the name of the default `Residue` in a `DisorderedResidue`. The default is
the first name read from the PDB file.
"""
defaultresname(disordered_res::DisorderedResidue) = disordered_res.default

"""
Access the default `Residue` in a `DisorderedResidue`. The default is the first
name read from the PDB file.
"""
defaultresidue(disordered_res::DisorderedResidue) = disordered_res.names[defaultresname(disordered_res)]

"""
Get the residue names in a `DisorderedResidue`. Default residue name is placed
first, then others are ordered alphabetically.
"""
resnames(disordered_res::DisorderedResidue) = sort(collect(keys(disordered_res.names)), lt= (res_name_one, res_name_two) -> (isless(res_name_one, res_name_two) && res_name_two != defaultresname(disordered_res)) || res_name_one == defaultresname(disordered_res))

resname(disordered_res::DisorderedResidue) = defaultresname(disordered_res)
resnumber(disordered_res::DisorderedResidue) = resnumber(defaultresidue(disordered_res))
inscode(disordered_res::DisorderedResidue) = inscode(defaultresidue(disordered_res))
ishetres(disordered_res::DisorderedResidue) = ishetres(defaultresidue(disordered_res))
atomnames(disordered_res::DisorderedResidue) = atomnames(defaultresidue(disordered_res))
atoms(disordered_res::DisorderedResidue) = atoms(defaultresidue(disordered_res))

chain(disordered_res::DisorderedResidue) = chain(defaultresidue(disordered_res))
model(disordered_res::DisorderedResidue) = model(defaultresidue(disordered_res))
structure(disordered_res::DisorderedResidue) = structure(defaultresidue(disordered_res))

isdisorderedres(::DisorderedResidue) = true

# Constructor acts as a setter for the default residue name
function DisorderedResidue(disordered_res::DisorderedResidue, default::AbstractString)
    @assert default in resnames(disordered_res) "The new default residue name must be present in the residue"
    return DisorderedResidue(disordered_res.names, default)
end


# Chain getters/setters
chainid(chain::Chain) = chain.id

"Get the sorted list of `AbstractResidue`s in a `Chain`."
resids(chain::Chain) = chain.res_list

"Access the dictionary of `AbstractResidue`s in a `Chain`."
residues(chain::Chain) = chain.residues

model(chain::Chain) = chain.model
structure(chain::Chain) = structure(model(chain))


# Model getters/setters

"Get the model number of a `Model`."
modelnumber(model::Model) = model.number

"Get the sorted chain IDs of the chains in a `Model` or `ProteinStructure`."
chainids(model::Model) = sortchainids(collect(keys(model.chains)))

"Access the dictionary of `Chain`s in a `Model`."
chains(model::Model) = model.chains

structure(model::Model) = model.structure


# ProteinStructure getters/setters

"Get the name of a `ProteinStructure`."
structurename(struc::ProteinStructure) = struc.name

"Get the sorted `Model` numbers from a `ProteinStructure`."
modelnumbers(struc::ProteinStructure) = sort(collect(keys(struc.models)))

"Access the dictionary of `Model`s in a `ProteinStructure`."
models(struc::ProteinStructure) = struc.models

"""
Get the default `Model` in a `ProteinStructure`. This is the `Model` with the
lowest model number.
"""
defaultmodel(struc::ProteinStructure) = struc.models[modelnumbers(struc)[1]]

chainids(struc::ProteinStructure) = countmodels(struc) > 0 ? chainids(defaultmodel(struc)) : Char[]


# Sorting functions

# Sort atoms by serial
Base.isless(atom_one::AbstractAtom, atom_two::AbstractAtom) = isless(serial(atom_one), serial(atom_two))

# Sort residues by chain, then hetero, then resumber, then ins code
function Base.isless(res_one::AbstractResidue, res_two::AbstractResidue)
    if chainidisless(chainid(res_one), chainid(res_two))
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


"""
Determine if the second residue follows the first in sequence. For this to be
true the residues need to have the same chain ID, both need to be
standard/hetero residues and the residue number of the second needs to be one
greater than that of the first.
"""
function sequentialresidues(res_first::AbstractResidue, res_second::AbstractResidue)
    if resnumber(res_second) == resnumber(res_first) + 1 && chainid(res_second) == chainid(res_first) && ishetres(res_second) == ishetres(res_first)
        return true
    else
        return false
    end
end


"""
Sort chain IDs. Chains are ordered by character sorting except the empty chain
ID comes last.
"""
sortchainids(chain_ids::Vector{Char}) = sort(chain_ids, lt=chainidisless)

"Determine if one chain ID should be sorted before another."
function chainidisless(chain_id_one::Char, chain_id_two::Char)
    if Int(chain_id_one) < Int(chain_id_two) && chain_id_one != ' ' && chain_id_two != ' '
        return true
    elseif chain_id_two == ' ' && chain_id_one != ' '
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
collectresidues(struc::ProteinStructure, selector_functions::Function...) = countmodels(struc) > 0 ? collectresidues(defaultmodel(struc), selector_functions...) : AbstractResidue[]
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
collectatoms(struc::ProteinStructure, selector_functions::Function...) = countmodels(struc) > 0 ? collectatoms(defaultmodel(struc), selector_functions...) : AbstractAtom[]
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
countchains(struc::ProteinStructure) = countmodels(struc) > 0 ? length(defaultmodel(struc)) : 0
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
`Model` and a `Vector{Model}` becomes a `ProteinStructure`.
"""
function organise(models::Vector{Model}; structure_name::AbstractString="")
    # Organise a Vector{Model} into a ProteinStructure
    struc = ProteinStructure(structure_name)
    for model in models
        @assert !(modelnumber(model) in modelnumbers(struc)) "Multiple models with the same model number found - cannot organise into a protein structure"
        struc[modelnumber(model)] = model
    end
    return struc
end

# Organise a Vector{Chain} into a Model
function organise(chains::Vector{Chain}; model_number::Integer=1)
    model = Model(model_number)
    for chain in chains
        @assert !(chainid(chain) in chainids(model)) "Multiple chains with the same chain ID found - cannot organise into a model"
        model[chainid(chain)] = chain
    end
    return model
end

# Organise a Vector{AbstractResidue} into a Vector{Chain}
function organise{T <: AbstractResidue}(residues::Vector{T})
    chains = Dict{Char, Dict{String, AbstractResidue}}()
    for res in residues
        chain_id = chainid(res)
        !haskey(chains, chain_id) ? chains[chain_id] = Dict{String, AbstractResidue}() : nothing
        res_id = resid(res)
        @assert !haskey(chains[chain_id], res_id) "Multiple residues with the same residue ID found - cannot organise into chains"
        chains[chain_id][res_id] = res
    end
    return [Chain(chain_id, map(resid, sort(collect(values(chains[chain_id])))), chains[chain_id]) for chain_id in sortchainids(collect(keys(chains)))]
end

# Organise a Vector{AbstractAtom} into a Vector{AbstractResidue}
function organise{T <: AbstractAtom}(atoms::Vector{T})
    # Key is chain ID, value is Dict where key is residue ID and value is list of atoms
    residues = Dict{Char, Dict{String, Vector{AbstractAtom}}}()
    # Key is chain ID, value is Dict where key is residue ID, value is Dict where key is residue name and value is list of atoms
    disordered_residues = Dict{Char, Dict{String, Dict{String, Vector{AbstractAtom}}}}()
    # Key is chain ID, value is Dict where key is residue ID and value is default residue name
    defaults = Dict{Char, Dict{String, String}}()
    for atom in atoms
        # Create chain Dict if required
        chain_id = chainid(atom)
        if !haskey(residues, chain_id)
            residues[chain_id] = Dict{String, Vector{AbstractAtom}}()
            disordered_residues[chain_id] = Dict{String, Dict{String, Vector{AbstractAtom}}}()
            defaults[chain_id] = Dict{String, String}()
        end
        res_id = resid(atom)
        # The residue ID is not yet present
        if !haskey(residues[chain_id], res_id) && !haskey(disordered_residues[chain_id], res_id)
            residues[chain_id][res_id] = AbstractAtom[]
            current_residue = residues[chain_id][res_id]
        # The disordered residue container could already exist
        elseif haskey(disordered_residues[chain_id], res_id)
            # If the res name isn't present, need to add a new Residue
            !haskey(disordered_residues[chain_id][res_id], resname(atom)) ? disordered_residues[chain_id][res_id][resname(atom)] = AbstractAtom[] : nothing
            current_residue = disordered_residues[chain_id][res_id][resname(atom)]
        # The disorered residue container doesn't exist and needs creating
        elseif resname(atom) != resname(residues[chain_id][res_id][1])
            disordered_residues[chain_id][res_id] = Dict(
                resname(residues[chain_id][res_id][1]) => residues[chain_id][res_id],
                resname(atom) => AbstractAtom[]
            )
            # The default res name is the first added
            defaults[chain_id][res_id] = resname(residues[chain_id][res_id][1])
            # Remove Residue now we have created a DisorderedResidue
            delete!(residues[chain_id], res_id)
            current_residue = disordered_residues[chain_id][res_id][resname(atom)]
        else
            current_residue = residues[chain_id][res_id]
        end
        atom_name = atomname(atom)
        index_found = 0
        for i in eachindex(current_residue)
            atom_name == atomname(current_residue[i]) ? (index_found = i; break) : nothing
        end
        @assert index_found == 0 "Multiple atoms with the same atom name on the same residue - cannot organise into residues:\n$(current_residue[index_found])\n$(atom)"
        push!(current_residue, atom)
    end
    residues_out = AbstractResidue[]
    for chain_id in sortchainids(collect(keys(residues)))
        residues_unord = AbstractResidue[]
        for res_id in keys(residues[chain_id])
            push!(residues_unord, Residue(residues[chain_id][res_id]))
        end
        for res_id in keys(disordered_residues[chain_id])
            new_disordered_res = DisorderedResidue(Dict(), defaults[chain_id][res_id])
            for res_name in keys(disordered_residues[chain_id][res_id])
                new_disordered_res.names[res_name] = Residue(disordered_residues[chain_id][res_id][res_name])
            end
            push!(residues_unord, new_disordered_res)
        end
        append!(residues_out, sort(residues_unord))
    end
    return residues_out
end

organise(model::Model; structure_name::AbstractString="") = organise([model]; structure_name=structure_name)
organise(chain::Chain; model_number::Integer=1) = organise([chain]; model_number=model_number)
organise(res::AbstractResidue) = organise(AbstractResidue[res])
organise(atom::AbstractAtom) = organise(AbstractAtom[atom])


"""
Form a list of `AbstractAtom`s from a list of `Atom`s. Combines disordered
atoms into disordered atom containers unless `remove_disorder` is `true`, in
which case removes all but the default location for disordered atoms.
"""
function formatomlist(atoms::Vector{Atom}; remove_disorder::Bool=false)
    # Key is (residue ID, residue name, atom name)
    atom_dic = Dict{Tuple{String, String, String}, AbstractAtom}()
    for atom in atoms
        atom_id = atomid(atom)
        # The atom does not exist so we can create it
        if !haskey(atom_dic, atom_id)
            atom_dic[atom_id] = atom
        # A disordered atom container already exists and the alt loc ID is not taken
        elseif isa(atom_dic[atom_id], DisorderedAtom) && !(altlocid(atom) in altlocids(atom_dic[atom_id]))
            # Add the new atom to the disordered atom container
            atom_dic[atom_id][altlocid(atom)] = atom
            # If the default alt loc requires changing, change it
            choosedefaultaltlocid(defaultatom(atom_dic[atom_id]), atom) != defaultaltlocid(atom_dic[atom_id]) ? atom_dic[atom_id] = DisorderedAtom(atom_dic[atom_id], altlocid(atom)) : nothing
        # The atom already exists and the alt loc IDs are different
        elseif isa(atom_dic[atom_id], Atom) && altlocid(atom) != altlocid(atom_dic[atom_id])
            # If we are removing disorder and the new atom is preferred to the old one, replace the old one
            if remove_disorder
                choosedefaultaltlocid(atom, atom_dic[atom_id]) == altlocid(atom) ? atom_dic[atom_id] = atom : nothing
            # If we are not removing disorder, create a new disordered atom container and add both atoms
            else
                atom_dic[atom_id] = DisorderedAtom(Dict(
                    altlocid(atom) => atom,
                    altlocid(atom_dic[atom_id]) => atom_dic[atom_id]
                ), choosedefaultaltlocid(atom, atom_dic[atom_id]))
            end
        else
            error("Two copies of the same atom have the same alternative location ID:\n$(atom_dic[atom_id])\n$(atom)\nformatomlist cannot be used on atom lists containing multiple models. In addition names are stripped of whitespace so identical names with different spacing, e.g. ' CA ' and 'CA  ', are not currently supported.")
        end
    end
    return sort(collect(values(atom_dic)))
end


# either make unsafe or deal with lists
#function addatomtomodel!(model::Model, serial::Int, name::Compat.ASCIIString, alt_loc_id::Char, coords::Vector{Float64}, occupancy::Float64, temp_fac::Float64, element::Compat.ASCIIString, charge::Compat.ASCIIString, het_atom::Bool, res_name::Compat.ASCIIString, chain_id::Char, res_number::Int, ins_code::Char; remove_disorder::Bool=false)
function addatomtomodel!(model::Model, line::Compat.ASCIIString, line_number::Integer=1; remove_disorder::Bool=false)
    het_atom = line[1] == 'H' # This is okay as the line has already been checked as an ATOM/HETATM record
    serial = parsestrict(line, 7, 11, Int, "Could not read atom serial number", line_number)
    atom_name = parsestrict(line, 13, 16, Compat.ASCIIString, "Could not read atom name", line_number) # Not stripped here for speed
    alt_loc_id = parsestrict(line, 17, 17, Char, "Could not read alt loc identifier", line_number)
    res_name = parsestrict(line, 18, 20, Compat.ASCIIString, "Could not read residue name", line_number) # Not stripped here for speed
    chain_id = parsestrict(line, 22, 22, Char, "Could not read chain ID", line_number)
    res_number = parsestrict(line, 23, 26, Int, "Could not read residue number", line_number)
    ins_code = parsestrict(line, 27, 27, Char, "Could not read insertion code", line_number)
    coords = [
        parsestrict(line, 31, 38, Float64, "Could not read x coordinate", line_number),
        parsestrict(line, 39, 46, Float64, "Could not read y coordinate", line_number),
        parsestrict(line, 47, 54, Float64, "Could not read z coordinate", line_number)
    ]
    occupancy = parselenient(line, 55, 60, Float64, 1.0)
    temp_fac = parselenient(line, 61, 66, Float64, 0.0)
    element = parselenient(line, 77, 78, Compat.ASCIIString, "") # Not stripped here for speed
    charge = parselenient(line, 79, 80, Compat.ASCIIString, "") # Not stripped here for speed

    # Add chain to model if necessary
    !haskey(chains(model), chain_id) ? model[chain_id] = Chain(chain_id, model) : nothing
    chain = model[chain_id]
    res_id = resid(het_atom, res_number, ins_code)
    # If residue does not exist in the chain, create a Residue
    if !haskey(residues(chain), res_id)
        chain[res_id] = Residue(res_name, res_number, ins_code, het_atom, chain)
        residue = chain[res_id]
    elseif isa(chain[res_id], Residue)
        # Residue exists in the chain and the residue names match
        # Add to that Residue
        if resname(chain[res_id]) == res_name
            residue = chain[res_id]
        # Residue exists in the chain but the residue names do not match
        # Create a DisorderedResidue
        else
            chain[res_id] = DisorderedResidue(Dict(
                resname(chain[res_id]) => chain[res_id],
                res_name => Residue(res_name, res_number, ins_code, het_atom, chain)
            ), resname(chain[res_id]))
            residue = disorderedres(chain[res_id], res_name)
        end
    else
        # DisorderedResidue exists in the chain and the residue names match
        # Add to that DisorderedResidue
        if res_name in resnames(chain[res_id])
            residue = disorderedres(chain[res_id], res_name)
        # DisorderedResidue exists in the chain and the residue names do not match
        # Create a new Residue in the DisorderedResidue
        else
            chain[res_id].names[res_name] = Residue(res_name, res_number, ins_code, het_atom, chain) # Make residues function for this?
            residue = disorderedres(chain[res_id], res_name)
        end
    end
    atom = Atom(serial, atom_name, alt_loc_id, coords, occupancy, temp_fac, element, charge, residue)
    # If atom does not exist in the residue, create an Atom
    if !haskey(atoms(residue), atom_name)
        residue[atom_name] = atom
    # Atom exists in the residue, atom names match and alt loc IDs are different
    elseif isa(residue[atom_name], Atom) && alt_loc_id != altlocid(residue[atom_name])
        # If we are removing disorder and the new atom is preferred to the old one, replace the old one
        if remove_disorder && choosedefaultaltlocid(atom, residue[atom_name]) == alt_loc_id
             residue[atom_name] = atom
        # If we are not removing disorder, create a new disordered atom container and add both atoms
        elseif !remove_disorder
            residue[atom_name] = DisorderedAtom(Dict(
                alt_loc_id => atom,
                altlocid(residue[atom_name]) => residue[atom_name]
            ), choosedefaultaltlocid(atom, residue[atom_name]))
        end
    # A disordered atom container already exists and the alt loc ID is not taken
    elseif isa(residue[atom_name], DisorderedAtom) && !(alt_loc_id in altlocids(residue[atom_name]))
        # Add the new atom to the disordered atom container
        residue[atom_name][alt_loc_id] = atom
        # If the default alt loc requires changing, change it
        if choosedefaultaltlocid(defaultatom(residue[atom_name]), atom) != defaultaltlocid(residue[atom_name])
            residue[atom_name] = DisorderedAtom(residue[atom_name], alt_loc_id)
        end
    else
        error("Two copies of the same atom have the same alternative location ID:\n$(residue[atom_name])\n$(atom)")
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
Organise elements into a `Model`. The keyword argument `model_number` sets the
model number (default `1`).
"""
organisemodel(chains::Vector{Chain}; model_number::Integer=1) = organise(chains; model_number=model_number)
organisemodel{T <: AbstractResidue}(residues::Vector{T}; model_number::Integer=1) = organise(organise(residues); model_number=model_number)
organisemodel{T <: AbstractAtom}(atoms::Vector{T}; model_number::Integer=1) = organise(organise(organise(atoms)); model_number=model_number)
organisemodel(chain::Chain; model_number::Integer=1) = organisemodel([chain]; model_number=model_number)
organisemodel(res::AbstractResidue; model_number::Integer=1) = organisemodel(AbstractResidue[res]; model_number=model_number)
organisemodel(atom::AbstractAtom; model_number::Integer=1) = organisemodel(AbstractAtom[atom]; model_number=model_number)


"""
Organise elements into a `ProteinStructure`. They keyword argument
`structure_name` sets the structure name (default `""`).
"""
organisestructure(models::Vector{Model}; structure_name::AbstractString="") = organise(models; structure_name=structure_name)
organisestructure(chains::Vector{Chain}; structure_name::AbstractString="", model_number::Integer=1) = organise(organisemodel(chains; model_number=model_number); structure_name=structure_name)
organisestructure{T <: AbstractResidue}(residues::Vector{T}; structure_name::AbstractString="", model_number::Integer=1) = organise(organisemodel(residues; model_number=model_number); structure_name=structure_name)
organisestructure{T <: AbstractAtom}(atoms::Vector{T}; structure_name::AbstractString="", model_number::Integer=1) = organise(organisemodel(atoms; model_number=model_number); structure_name=structure_name)
organisestructure(model::Model; structure_name::AbstractString="") = organisestructure([model]; structure_name=structure_name)
organisestructure(chain::Chain; structure_name::AbstractString="", model_number::Integer=1) = organisestructure([chain]; structure_name=structure_name, model_number=model_number)
organisestructure(res::AbstractResidue; structure_name::AbstractString="", model_number::Integer=1) = organisestructure(AbstractResidue[res]; structure_name=structure_name, model_number=model_number)
organisestructure(atom::AbstractAtom; structure_name::AbstractString="", model_number::Integer=1) = organisestructure(AbstractAtom[atom]; structure_name=structure_name, model_number=model_number)


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
