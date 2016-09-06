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
    AtomRecord,
    StructuralElementOrList,
    serial,
    atomname,
    altlocid,
    x,
    x!,
    y,
    y!,
    z,
    z!,
    coords,
    coords!,
    occupancy,
    tempfac,
    element,
    charge,
    residue,
    ishetatom,
    isdisorderedatom,
    defaultaltlocid,
    defaultatom,
    altlocids,
    atomid,
    resname,
    resnumber,
    inscode,
    ishetres,
    ishetero,
    resid,
    atomnames,
    atoms,
    isdisorderedres,
    disorderedres,
    defaultresname,
    defaultresidue,
    resnames,
    chain,
    chainid,
    resids,
    residues,
    model,
    modelnumber,
    chainids,
    chains,
    structure,
    structurename,
    modelnumbers,
    models,
    defaultmodel,
    sequentialresidues,
    applyselectors,
    applyselectors!,
    collectmodels,
    countmodels,
    collectchains,
    countchains,
    collectresidues,
    countresidues,
    collectatoms,
    countatoms,
    choosedefaultaltlocid,
    stdatomselector,
    hetatomselector,
    atomnameselector,
    calphaatomnames,
    calphaselector,
    cbetaatomnames,
    cbetaselector,
    backboneatomnames,
    backboneselector,
    heavyatomselector,
    resnameselector,
    waterresnames,
    waterselector,
    stdresselector,
    hetresselector,
    disorderselector,
    hydrogenselector,
    AminoAcidSequence


using Bio.Seq
import Bio.Seq.AminoAcidSequence
import Bio.Seq.threeletter_to_aa


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

"A container to hold different locations of the same atom."
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
file.
"""
immutable AtomRecord
    het_atom::Bool
    serial::Int
    atom_name::String
    alt_loc_id::Char
    res_name::String
    chain_id::Char
    res_number::Int
    ins_code::Char
    coords::Vector{Float64}
    occupancy::Float64
    temp_fac::Float64
    element::String
    charge::String
end


# Constructors without sub-elements and without super-elements

ProteinStructure(name::AbstractString) = ProteinStructure(name, Dict())

ProteinStructure() = ProteinStructure("")


Model(number::Integer, struc::ProteinStructure) = Model(number, Dict(), struc)

Model(struc::ProteinStructure) = Model(1, struc)

Model(number::Integer) = Model(number, ProteinStructure())

Model() = Model(ProteinStructure())


Chain(id::Char, mod::Model) = Chain(id, [], Dict(), mod)

Chain(mod::Model) = Chain('A', mod)

Chain(id::Char) = Chain(id, Model())

Chain() = Chain('A')


function Residue(name::AbstractString,
                number::Integer,
                ins_code::Char,
                het_res::Bool,
                ch::Chain)
    return Residue(name, number, ins_code, het_res, [], Dict(), ch)
end

Residue(ch::Chain) = Residue("ALA", 1, ' ', false, ch)

function Residue(name::AbstractString,
                number::Integer,
                ins_code::Char,
                het_res::Bool)
    return Residue(name, number, ins_code, het_res, Chain())
end

Residue() = Residue("ALA", 1, ' ', false)


function Atom(serial::Integer,
                name::AbstractString,
                alt_loc_id::Char,
                coords::Vector{Float64},
                occupancy::Float64,
                temp_fac::Float64,
                element::AbstractString,
                charge::AbstractString)
    return Atom(serial, name, alt_loc_id, coords, occupancy, temp_fac, element, charge, Residue())
end

Atom(res::Residue) = Atom(1, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "  ", "  ", res)

Atom() = Atom(1, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "  ", "  ")


"A `StructuralElement` or `Vector` of `StructuralElement`s."
typealias StructuralElementOrList Union{
        StructuralElement,
        Vector{Model},
        Vector{Chain},
        Vector{AbstractResidue},
        Vector{Residue},
        Vector{DisorderedResidue},
        Vector{AbstractAtom},
        Vector{Atom},
        Vector{DisorderedAtom}
    }


# Allow accessing sub elements contained in an element like a dictionary
# e.g. allows you to do res[atom_name] rather than res.atoms[atom_name]
# setindex! should be used with caution as it can lead to inconsistencies
# e.g. adding an atom to a residue atom dict does not update the atom list
# Hence it is used internally but not documented for general use

# Accessing a DisorderedAtom with a character returns the Atom with that alt
#   loc ID
Base.getindex(dis_at::DisorderedAtom, alt_loc_id::Char) = dis_at.alt_loc_ids[alt_loc_id]

function Base.setindex!(dis_at::DisorderedAtom, at::Atom, alt_loc_id::Char)
    dis_at.alt_loc_ids[alt_loc_id] = at
    return dis_at
end

# Accessing a Residue with an AbstractString returns the AbstractAtom with that
#   atom name
Base.getindex(res::Residue, atom_name::AbstractString) = findatombyname(res, atom_name)

function Base.setindex!(res::Residue, at::AbstractAtom, atom_name::AbstractString)
    res.atoms[atom_name] = at
    return res
end

# Accessing a DisorderedResidue with an AbstractString returns the AbstractAtom
#   in the default Residue with that atom name
# This is not necessarily intuitive, it may be expected to return the Residue
#   with that residue name
# However this way accessing an AbstractResidue always returns an AbstractAtom
Base.getindex(dis_res::DisorderedResidue, atom_name::AbstractString) = findatombyname(defaultresidue(dis_res), atom_name)

function Base.setindex!(dis_res::DisorderedResidue, at::AbstractAtom, atom_name::AbstractString)
    defaultresidue(dis_res)[atom_name] = at
    return dis_res
end

# Accessing a Chain with an AbstractString returns the AbstractResidue with that
#  residue ID
Base.getindex(ch::Chain, res_id::AbstractString) = ch.residues[res_id]

function Base.setindex!(ch::Chain, res::AbstractResidue, res_id::AbstractString)
    ch.residues[res_id] = res
    return ch
end

# Accessing a Chain with an Integer returns the AbstractResidue with that residue ID
#   converted to a string
Base.getindex(ch::Chain, res_n::Integer) = ch.residues[string(res_n)]

function Base.setindex!(ch::Chain, res::AbstractResidue, res_n::Integer)
    ch.residues[string(res_n)] = res
    return ch
end

# Accessing a Model with a Char returns the Chain with that chain ID
Base.getindex(mod::Model, ch_id::Char) = mod.chains[ch_id]

function Base.setindex!(mod::Model, ch::Chain, ch_id::Char)
    mod.chains[ch_id] = ch
    return mod
end

# Accessing a ProteinStructure with an Integer returns the Model with that model
#   number
Base.getindex(struc::ProteinStructure, mod_n::Integer) = struc.models[mod_n]

function Base.setindex!(struc::ProteinStructure, mod::Model, mod_n::Integer)
    struc.models[mod_n] = mod
    return struc
end

# Accessing a ProteinStructure with a Char returns the Chain with that chain ID
#   on the default model
Base.getindex(struc::ProteinStructure, ch_id::Char) = defaultmodel(struc)[ch_id]

function Base.setindex!(struc::ProteinStructure, ch::Chain, ch_id::Char)
    defaultmodel(struc)[ch_id] = ch
    return struc
end

# Check if an atom name exists in a residue as a whitespace-padded version
function findatombyname(res::Residue, atom_name::AbstractString)
    # Look for atom name directly
    if atom_name in keys(res.atoms)
        return res.atoms[atom_name]
    # Pad out name to 4 characters to read PDB atom names with whitespace
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
    # Could not find atom name
    throw(KeyError(atom_name))
end


# Getters and setters for structural elements

"Get the serial number of an `AbstractAtom`."
serial(at::Atom) = at.serial
serial(dis_at::DisorderedAtom) = serial(defaultatom(dis_at))


"""
Get the atom name of an `AbstractAtom`.
`spaces` determines whether surrounding whitespace is retained (default
`false`).
"""
function atomname(at::Atom; spaces::Bool=false)
    if spaces
        return at.name
    else
        return strip(at.name)
    end
end

atomname(dis_at::DisorderedAtom; spaces::Bool=false) = atomname(defaultatom(dis_at), spaces=spaces)


"Get the alternative location ID of an `AbstractAtom`."
altlocid(at::Atom) = at.alt_loc_id
altlocid(dis_at::DisorderedAtom) = defaultaltlocid(dis_at)


"Get the x coordinate of an `AbstractAtom`."
x(at::Atom) = at.coords[1]
x(dis_at::DisorderedAtom) = x(defaultatom(dis_at))

"""
Set the x coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
x!(at::Atom, x::Real) = (at.coords[1] = x; at)
x!(dis_at::DisorderedAtom, x::Real) = x!(defaultatom(dis_at), x)


"Get the y coordinate of an `AbstractAtom`."
y(at::Atom) = at.coords[2]
y(dis_at::DisorderedAtom) = y(defaultatom(dis_at))

"""
Set the y coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
y!(at::Atom, y::Real) = (at.coords[2] = y; at)
y!(dis_at::DisorderedAtom, y::Real) = y!(defaultatom(dis_at), y)


"Get the z coordinate of an `AbstractAtom`."
z(at::Atom) = at.coords[3]
z(dis_at::DisorderedAtom) = z(defaultatom(dis_at))

"""
Set the z coordinate of an `AbstractAtom`. For `DisorderedAtom`s only the
default atom is updated.
"""
z!(at::Atom, z::Real) = (at.coords[3] = z; at)
z!(dis_at::DisorderedAtom, z::Real) = z!(defaultatom(dis_at), z)


"Get the atomic coordinates of an `AbstractAtom` as a `Vector{Float64}`."
coords(at::Atom) = at.coords
coords(dis_at::DisorderedAtom) = coords(defaultatom(dis_at))

"""
Set the coordinates of an `AbstractAtom` to a `Vector{Float64}` of length 3. For
`DisorderedAtom`s only the default atom is updated.
"""
function coords!(at::Atom, new_coords::Vector{Float64})
    @assert length(new_coords) == 3 "3 coordinates must be given"
    x!(at, new_coords[1])
    y!(at, new_coords[2])
    z!(at, new_coords[3])
    return at
end

function coords!(dis_at::DisorderedAtom, coords::Vector{Float64})
    coords!(defaultatom(dis_at), coords)
end


"Get the occupancy of an `AbstractAtom`. Defaults to `1.0`."
occupancy(at::Atom) = at.occupancy
occupancy(dis_at::DisorderedAtom) = occupancy(defaultatom(dis_at))


"Get the temperature factor of an `AbstractAtom`. Defaults to `0.0`."
tempfac(at::Atom) = at.temp_fac
tempfac(dis_at::DisorderedAtom) = tempfac(defaultatom(dis_at))


"""
Get the element of an `AbstractAtom`. Defaults to `\"  \"`.
`spaces` determines whether surrounding whitespace is retained (default
`false`).
"""
function element(at::Atom; spaces::Bool=false)
    if spaces
        return at.element
    else
        return strip(at.element)
    end
end

element(dis_at::DisorderedAtom; spaces::Bool=false) = element(defaultatom(dis_at), spaces=spaces)


"""
Get the charge on an `AbstractAtom`. Defaults to `\"  \"`.
`spaces` determines whether surrounding whitespace is retained (default
`false`).
"""
function charge(at::Atom; spaces::Bool=false)
    if spaces
        return at.charge
    else
        return strip(at.charge)
    end
end

charge(dis_at::DisorderedAtom; spaces::Bool=false) = charge(defaultatom(dis_at), spaces=spaces)


"Get the `Residue` that an `AbstractAtom` belongs to."
residue(at::Atom) = at.residue
residue(dis_at::DisorderedAtom) = residue(defaultatom(dis_at))
residue(res::AbstractResidue) = res


"""
Determines if an `AbstractAtom` represents a non-hetero atom, e.g. came from an
ATOM record in a Protein Data Bank (PDB) file.
"""
ishetatom(at::Atom) = ishetres(residue(at))
ishetatom(dis_at::DisorderedAtom) = ishetatom(defaultatom(dis_at))


"""
Determines if an `AbstractAtom` is a `DisorderedAtom`, i.e. if there are
multiple locations present for an atom.
"""
isdisorderedatom(at::AbstractAtom) = isa(at, DisorderedAtom)


"""
Get the alternative location ID of the default `Atom` in a `DisorderedAtom`.
The default is the highest occupancy, or lowest character alternative location
ID for ties.
"""
defaultaltlocid(dis_at::DisorderedAtom) = dis_at.default


# Constructor acts as a setter for the default alt loc ID
function DisorderedAtom(dis_at::DisorderedAtom, default::Char)
    @assert default in altlocids(dis_at) "The new default alternative location ID must be present in the atom"
    return DisorderedAtom(dis_at.alt_loc_ids, default)
end


"""
Access the default `Atom` in a `DisorderedAtom`.
The default is the highest occupancy, or lowest character alternative location
ID for ties.
"""
defaultatom(dis_at::DisorderedAtom) = dis_at[defaultaltlocid(dis_at)]


"""
Get the list of alternative location IDs in an `AbstractAtom`, sorted by atom
serial.
"""
function altlocids(dis_at::DisorderedAtom)
    return sort(collect(keys(dis_at.alt_loc_ids)),
        by= alt_loc_id -> dis_at[alt_loc_id])
end

altlocids(at::Atom) = [altlocid(at)]


"""
Get a descriptive atom ID for an `AbstractAtom` as a `Tuple` of the form
(full residue ID, residue name, atom name).
"""
atomid(at::Atom) = (resid(at, full=true), resname(at), atomname(at))
atomid(dis_at::DisorderedAtom) = atomid(defaultatom(dis_at))


"""
Get the residue name of an `AbstractAtom` or `AbstractResidue`.
`spaces` determines whether surrounding whitespace is retained (default
`false`).
"""
function resname(res::Residue; spaces::Bool=false)
    if spaces
        return res.name
    else
        return strip(res.name)
    end
end

resname(at::Atom; spaces::Bool=false) = resname(residue(at), spaces=spaces)
resname(dis_at::DisorderedAtom; spaces::Bool=false) = resname(defaultatom(dis_at), spaces=spaces)
resname(dis_res::DisorderedResidue; spaces::Bool=false) = resname(defaultresidue(dis_res), spaces=spaces)


"Get the residue number of an `AbstractAtom` or `AbstractResidue`."
resnumber(at::Atom) = resnumber(residue(at))
resnumber(dis_at::DisorderedAtom) = resnumber(defaultatom(dis_at))
resnumber(res::Residue) = res.number
resnumber(dis_res::DisorderedResidue) = resnumber(defaultresidue(dis_res))


"Get the insertion code of an `AbstractAtom` or `AbstractResidue`."
inscode(at::Atom) = inscode(residue(at))
inscode(dis_at::DisorderedAtom) = inscode(defaultatom(dis_at))
inscode(res::Residue) = res.ins_code
inscode(dis_res::DisorderedResidue) = inscode(defaultresidue(dis_res))


"""
Determines if an `AbstractResidue` represents a hetero molecule, e.g. consists
of HETATM records from a Protein Data Bank (PDB) file.
"""
ishetres(res::Residue) = res.het_res
ishetres(dis_res::DisorderedResidue) = ishetres(defaultresidue(dis_res))


"""
Determines if an `AbstractAtom` represents a hetero atom, e.g. came from a
HETATM record in a Protein Data Bank (PDB) file, or if an `AbstractResidue`
represents a hetero molecule, e.g. consists of HETATM records from a PDB file.
"""
ishetero(at::AbstractAtom) = ishetatom(at)
ishetero(res::AbstractResidue) = ishetres(res)


"""
Get a descriptive residue ID string for an `AbstractAtom` or `AbstractResidue`.
Format is residue number then insertion code with \"H_\" in front for hetero
residues. If `full` equals `true` the chain ID is also added after a colon
(default `false`). Examples are \"50A\", \"H_20\" and \"10:A\".
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

resid(at::AbstractAtom; full::Bool=false) = resid(residue(at), full=full)

# Shortcut from atomic properties
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


"""
Get the sorted list of `AbstractAtom`s in an `AbstractResidue`.
`spaces` determines whether surrounding whitespace is retained (default
`false`).
"""
function atomnames(res::Residue; spaces::Bool=false)
    if spaces
        return res.atom_list
    else
        return map(strip, res.atom_list)
    end
end

atomnames(dis_res::DisorderedResidue; spaces::Bool=false) = atomnames(defaultresidue(dis_res), spaces=spaces)


"Access the dictionary of `AbstractAtom`s in an `AbstractResidue`."
atoms(res::Residue) = res.atoms
atoms(dis_res::DisorderedResidue) = atoms(defaultresidue(dis_res))


"""
Determine if an `AbstractResidue` is a `DisorderedResidue`, i.e. there are
multiple residue names with the same residue ID.
"""
isdisorderedres(res::AbstractResidue) = isa(res, DisorderedResidue)


"Access the `Residue` in a `DisorderedResidue` with a given residue name."
disorderedres(dis_res::DisorderedResidue, res_name::AbstractString) = dis_res.names[res_name]


"""
Get the name of the default `Residue` in a `DisorderedResidue`. The default is
the first name read in.
"""
defaultresname(dis_res::DisorderedResidue) = dis_res.default


"""
Access the default `Residue` in a `DisorderedResidue`. The default is the first
name read in.
"""
defaultresidue(dis_res::DisorderedResidue) = dis_res.names[defaultresname(dis_res)]


"""
Get the residue names in an `AbstractResidue`. For a `DisorderedResidue` there
will be multiple residue names - in this case the default residue name is placed
first, then the others are ordered alphabetically.
"""
function resnames(dis_res::DisorderedResidue)
    return sort(collect(keys(dis_res.names)),
        lt= (res_name_one, res_name_two) ->
            (isless(res_name_one, res_name_two) && res_name_two != defaultresname(dis_res)) ||
            res_name_one == defaultresname(dis_res)
    )
end

resnames(res::Residue) = [resname(res, spaces=true)]


# Constructor acts as a setter for the default residue name
function DisorderedResidue(dis_res::DisorderedResidue, default::AbstractString)
    @assert default in resnames(dis_res) "The new default residue name must be present in the residue"
    return DisorderedResidue(dis_res.names, default)
end


"Get the `Chain` that an `AbstractAtom` or `AbstractResidue` belongs to."
chain(at::Atom) = chain(residue(at))
chain(dis_at::DisorderedAtom) = chain(defaultatom(dis_at))
chain(res::Residue) = res.chain
chain(dis_res::DisorderedResidue) = chain(defaultresidue(dis_res))
chain(ch::Chain) = ch


"Get the chain ID of an `AbstractAtom`, `AbstractResidue` or `Chain`."
chainid(el::Union{AbstractResidue, AbstractAtom}) = chainid(chain(el))
chainid(ch::Chain) = ch.id


"Get the sorted list of `AbstractResidue`s in a `Chain`."
resids(ch::Chain) = ch.res_list


"Access the dictionary of `AbstractResidue`s in a `Chain`."
residues(ch::Chain) = ch.residues


"""
Get the `Model` that an `AbstractAtom`, `AbstractResidue` or `Chain` belongs to.
"""
model(at::Atom) = model(chain(at))
model(dis_at::DisorderedAtom) = model(defaultatom(dis_at))
model(res::Residue) = model(chain(res))
model(dis_res::DisorderedResidue) = model(defaultresidue(dis_res))
model(ch::Chain) = ch.model
model(mod::Model) = mod


"""
Get the model number of a `Model`, `Chain`, `AbstractResidue` or `AbstractAtom`.
"""
modelnumber(mod::Model) = mod.number
modelnumber(el::Union{Chain, AbstractResidue, AbstractAtom}) = modelnumber(model(el))


"""
Get the sorted chain IDs of the chains in a `Model`, or the default `Model` of a
`ProteinStructure`.
"""
chainids(mod::Model) = map(chainid, sort(collect(values(chains(mod)))))

function chainids(struc::ProteinStructure)
    if countmodels(struc) > 0
        return chainids(defaultmodel(struc))
    else
        return Char[]
    end
end


"""
Access the dictionary of `Chain`s in a `Model`, or the default `Model` of a
`ProteinStructure`.
"""
chains(mod::Model) = mod.chains

function chains(struc::ProteinStructure)
    if countmodels(struc) > 0
        return chains(defaultmodel(struc))
    else
        return Dict{Char, Chain}()
    end
end


"""
Get the `ProteinStructure` that an `AbstractAtom`, `AbstractResidue`, `Chain`
or `Model` belongs to.
"""
structure(at::Atom) = structure(model(at))
structure(dis_at::DisorderedAtom) = structure(defaultatom(dis_at))
structure(res::Residue) = structure(model(res))
structure(dis_res::DisorderedResidue) = structure(defaultresidue(dis_res))
structure(ch::Chain) = structure(model(ch))
structure(mod::Model) = mod.structure
structure(struc::ProteinStructure) = struc


"Get the structure name of a `StructuralElement`."
structurename(el::Union{Model, Chain, AbstractResidue, AbstractAtom}) = structurename(structure(el))
structurename(struc::ProteinStructure) = struc.name


"Get the sorted model numbers from a `ProteinStructure`."
function modelnumbers(struc::ProteinStructure)
    return map(modelnumber, sort(collect(values(models(struc)))))
end


"Access the dictionary of `Model`s in a `ProteinStructure`."
models(struc::ProteinStructure) = struc.models


"""
Get the default `Model` in a `ProteinStructure`. This is the `Model` with the
lowest model number.
"""
defaultmodel(struc::ProteinStructure) = struc.models[modelnumbers(struc)[1]]


# Sort lists of elements

# Sort atoms by serial
function Base.isless(at_one::AbstractAtom, at_two::AbstractAtom)
    return isless(serial(at_one), serial(at_two))
end

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
function Base.isless(ch_one::Chain, ch_two::Chain)
    if Int(chainid(ch_one)) < Int(chainid(ch_two)) && chainid(ch_one) != ' ' && chainid(ch_two) != ' '
        return true
    elseif chainid(ch_two) == ' ' && chainid(ch_one) != ' '
        return true
    else
        return false
    end
end

# Sort models by model number
function Base.isless(mod_one::Model, mod_two::Model)
    return isless(modelnumber(mod_one), modelnumber(mod_two))
end


"""
Determine if the second residue follows the first in sequence. For this to be
`true` the residues need to have the same chain ID, both need to be
standard/hetero residues and the residue number of the second needs to be one
greater than that of the first (or the residue numbers the same and the
insertion code of the second greater than the first).
"""
function sequentialresidues(res_first::AbstractResidue, res_second::AbstractResidue)
    if chainid(res_second) == chainid(res_first) &&
            ishetres(res_second) == ishetres(res_first)
        if resnumber(res_second) == resnumber(res_first) + 1
            return true
        elseif resnumber(res_second) == resnumber(res_first) &&
                inscode(res_second) > inscode(res_first)
            return true
        end
    end
    return false
end


# Iterators to yield sub elements when looping over an element

# Iterating over a ProteinStructure yields Models
Base.length(struc::ProteinStructure) = length(modelnumbers(struc))
Base.start(::ProteinStructure) = 1
Base.next(struc::ProteinStructure, state) = (struc[modelnumbers(struc)[state]], state + 1)
Base.done(struc::ProteinStructure, state) = state > length(struc)
Base.eltype(::Type{ProteinStructure}) = Model

# Iterating over a Model yields Chains
Base.length(mod::Model) = length(chainids(mod))
Base.start(::Model) = 1
Base.next(mod::Model, state) = (mod[chainids(mod)[state]], state + 1)
Base.done(mod::Model, state) = state > length(mod)
Base.eltype(::Type{Model}) = Chain

# Iterating over a Chain yields AbstractResidues
Base.length(ch::Chain) = length(resids(ch))
Base.start(::Chain) = 1
Base.next(ch::Chain, state) = (ch[resids(ch)[state]], state + 1)
Base.done(ch::Chain, state) = state > length(ch)
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
Base.length(dis_res::DisorderedResidue) = length(atomnames(dis_res))
Base.start(::DisorderedResidue) = 1
Base.next(dis_res::DisorderedResidue, state) = (defaultresidue(dis_res)[atomnames(dis_res)[state]], state + 1)
Base.done(dis_res::DisorderedResidue, state) = state > length(dis_res)
Base.eltype(::Type{DisorderedResidue}) = AbstractAtom

# Iterating over an Atom returns itself
# This is not necessarily intuitive, it may be expected to not be an iterator
# However this way iterating over an AbstractAtom always yields Atoms
Base.length(::Atom) = 1
Base.start(::Atom) = 1
Base.next(at::Atom, state) = (at, state + 1)
Base.done(::Atom, state) = state > 1
Base.eltype(::Type{Atom}) = Atom

# Iterating over a DisorderedAtom yields Atoms
Base.length(dis_at::DisorderedAtom) = length(altlocids(dis_at))
Base.start(::DisorderedAtom) = 1
Base.next(dis_at::DisorderedAtom, state) = (dis_at[altlocids(dis_at)[state]], state + 1)
Base.done(dis_at::DisorderedAtom, state) = state > length(dis_at)
Base.eltype(::Type{DisorderedAtom}) = Atom


"""
Returns a copy of a `Vector` of `StructuralElement`s with all elements that do
not return `true` from the selector functions removed.
"""
function applyselectors{T <: StructuralElement}(els::Vector{T}, selectors::Function...)
    new_list = copy(els)
    applyselectors!(new_list, selectors...)
    return new_list
end

"Runs `applyselectors` in place."
function applyselectors!{T <: StructuralElement}(els::Vector{T}, selectors::Function...)
    for selector in selectors
        filter!(selector, els)
    end
    return els
end


"""
Returns a sorted `Vector` of the models in a `StructuralElementOrList`.
Additional arguments are model selector functions - only models that return
`true` from the functions are retained.
"""
collectmodels(struc::ProteinStructure) = collect(struc)

collectmodels(mod::Model) = [mod]

collectmodels(el::Union{Chain, AbstractResidue, AbstractAtom}) = [model(el)]

collectmodels(mods::Vector{Model}) = sort(mods)

function collectmodels{T <: Union{Chain, AbstractResidue, AbstractAtom}}(els::Vector{T})
    mod_list = Model[]
    for el in els
        if !(model(el) in mod_list)
            push!(mod_list, model(el))
        end
    end
    return sort(mod_list)
end

# One selector explicitly defined to prevent this being called without selectors
function collectmodels(el::StructuralElementOrList,
                    model_selector::Function,
                    model_selectors::Function...)
    return applyselectors(collectmodels(el), model_selector, model_selectors...)
end


"""
Return the number of `Model`s in a `StructuralElementOrList`.
Additional arguments are model selector functions - only models that return
`true` from the functions are counted.
"""
function countmodels(el::StructuralElementOrList, model_selectors::Function...)
    return length(collectmodels(el, model_selectors...))
end

countmodels(struc::ProteinStructure) = length(struc)

countmodels(mods::Vector{Model}) = length(mods)


"""
Returns a sorted `Vector` of the chains in a `StructuralElementOrList`.
Additional arguments are chain selector functions - only chains that return
`true` from the functions are retained.
"""
function collectchains(struc::ProteinStructure)
    if countmodels(struc) > 0
        return collectchains(defaultmodel(struc))
    else
        return Chain[]
    end
end

collectchains(mod::Model) = collect(mod)

collectchains(ch::Chain) = [ch]

collectchains(el::Union{AbstractResidue, AbstractAtom}) = [chain(el)]

function collectchains(mods::Vector{Model})
    ch_list = Chain[]
    for mod in mods
        append!(ch_list, collectchains(mod))
    end
    return sort(ch_list)
end

collectchains(chs::Vector{Chain}) = sort(chs)

function collectchains{T <: Union{AbstractResidue, AbstractAtom}}(els::Vector{T})
    ch_list = Chain[]
    for el in els
        if !(chain(el) in ch_list)
            push!(ch_list, chain(el))
        end
    end
    return sort(ch_list)
end

function collectchains(el::StructuralElementOrList,
                    chain_selector::Function,
                    chain_selectors::Function...)
    return applyselectors(collectchains(el), chain_selector, chain_selectors...)
end


"""
Return the number of `Chain`s in a `StructuralElementOrList`.
Additional arguments are chain selector functions - only chains that return
`true` from the functions are counted.
"""
function countchains(el::StructuralElementOrList, chain_selectors::Function...)
    return length(collectchains(el, chain_selectors...))
end

countchains(mod::Model) = length(mod)

countchains(chs::Vector{Chain}) = length(chs)


"""
Returns a sorted `Vector` of the residues in a `StructuralElementOrList`.
Additional arguments are residue selector functions - only residues that return
`true` from the functions are retained.
"""
function collectresidues(struc::ProteinStructure)
    if countmodels(struc) > 0
        return collectresidues(defaultmodel(struc))
    else
        return AbstractResidue[]
    end
end

function collectresidues(el::Union{Model, Vector{Model}, Vector{Chain}})
    res_list = AbstractResidue[]
    for sub_el in el
        append!(res_list, collectresidues(sub_el))
    end
    return sort(res_list)
end

collectresidues(ch::Chain) = collect(ch)

collectresidues(res::AbstractResidue) = AbstractResidue[res]

collectresidues(at::AbstractAtom) = AbstractResidue[residue(at)]

# Note output is always Vector{AbstractResidue} unless input was Vector{Residue}
# or Vector{DisorderedResidue}, in which case output is same type as input type
collectresidues{T <: AbstractResidue}(res_list::Vector{T}) = sort(res_list)

function collectresidues{T <: AbstractAtom}(at_list::Vector{T})
    res_list = AbstractResidue[]
    for at in at_list
        if !(residue(at) in res_list)
            push!(res_list, residue(at))
        end
    end
    return sort(res_list)
end

function collectresidues(el::StructuralElementOrList,
                    residue_selector::Function,
                    residue_selectors::Function...)
    return applyselectors(collectresidues(el), residue_selector, residue_selectors...)
end


"""
Return the number of residues in a `StructuralElementOrList`.
Additional arguments are residue selector functions - only residues that return
`true` from the functions are counted.
"""
function countresidues(el::StructuralElementOrList, residue_selectors::Function...)
    return length(collectresidues(el, residue_selectors...))
end

countresidues(ch::Chain) = length(ch)

function countresidues(res_list::Union{Vector{AbstractResidue},
        Vector{Residue}, Vector{DisorderedResidue}})
    return length(res_list)
end


"""
Returns a sorted `Vector` of the atoms in a `StructuralElementOrList`.
Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function collectatoms(struc::ProteinStructure)
    if countmodels(struc) > 0
        return collectatoms(defaultmodel(struc))
    else
        return AbstractAtom[]
    end
end

function collectatoms(el::Union{Model, Chain, Vector{Model}, Vector{Chain},
        Vector{AbstractResidue}, Vector{Residue}, Vector{DisorderedResidue}})
    at_list = AbstractAtom[]
    for sub_el in el
        append!(at_list, collectatoms(sub_el))
    end
    return sort(at_list)
end

collectatoms(res::AbstractResidue) = collect(res)

collectatoms(at::AbstractAtom) = AbstractAtom[at]

# Note output is always Vector{AbstractAtom} unless input was Vector{Atom} or
# Vector{DisorderedAtom}, in which case output is same type as input type
collectatoms{T <: AbstractAtom}(at_list::Vector{T}) = sort(at_list)

function collectatoms(el::StructuralElementOrList, atom_selector::Function, atom_selectors::Function...)
    return applyselectors(collectatoms(el), atom_selector, atom_selectors...)
end


"""
Return the number of atoms in a `StructuralElementOrList`.
Additional arguments are atom selector functions - only atoms that return
`true` from the functions are counted.
"""
function countatoms(el::StructuralElementOrList, atom_selectors::Function...)
    return length(collectatoms(el, atom_selectors...))
end

countatoms(res::AbstractResidue) = length(res)

function countatoms(at_list::Union{Vector{AbstractAtom},
        Vector{Atom}, Vector{DisorderedAtom}})
    return length(at_list)
end


# Add an atom represented in an AtomRecord to a Model
# Unsafe as sub-element lists are not updated (for speed)
# fixlists! should be run after all additions to update the sub-element lists
function unsafe_addatomtomodel!(mod::Model,
                    atom_rec::AtomRecord;
                    remove_disorder::Bool=false)
    # Add chain to model if necessary
    if !haskey(chains(mod), atom_rec.chain_id)
        mod[atom_rec.chain_id] = Chain(atom_rec.chain_id, mod)
    end
    ch = mod[atom_rec.chain_id]
    res_id = resid(atom_rec.het_atom, atom_rec.res_number, atom_rec.ins_code)
    # If residue does not exist in the chain, create a Residue
    if !haskey(residues(ch), res_id)
        ch[res_id] = Residue(
                    atom_rec.res_name,
                    atom_rec.res_number,
                    atom_rec.ins_code,
                    atom_rec.het_atom,
                    ch)
        res = ch[res_id]
    elseif isa(ch[res_id], Residue)
        # Residue exists in the chain and the residue names match
        # Add to that Residue
        if resname(ch[res_id]) == atom_rec.res_name
            res = ch[res_id]
        # Residue exists in the chain but the residue names do not match
        # Create a DisorderedResidue
        else
            ch[res_id] = DisorderedResidue(Dict(
                resname(ch[res_id]) => ch[res_id],
                atom_rec.res_name => Residue(
                    atom_rec.res_name,
                    atom_rec.res_number,
                    atom_rec.ins_code,
                    atom_rec.het_atom,
                    ch)
            ), resname(ch[res_id]))
            res = disorderedres(ch[res_id], atom_rec.res_name)
        end
    else
        # DisorderedResidue exists in the chain and the residue names match
        # Add to that DisorderedResidue
        if atom_rec.res_name in resnames(ch[res_id])
            res = disorderedres(ch[res_id], atom_rec.res_name)
        # DisorderedResidue exists in the chain and the residue names do not match
        # Create a new Residue in the DisorderedResidue
        else
            ch[res_id].names[atom_rec.res_name] = Residue(
                    atom_rec.res_name,
                    atom_rec.res_number,
                    atom_rec.ins_code,
                    atom_rec.het_atom,
                    ch)
            res = disorderedres(ch[res_id], atom_rec.res_name)
        end
    end
    at = Atom(
        atom_rec.serial,
        atom_rec.atom_name,
        atom_rec.alt_loc_id,
        atom_rec.coords,
        atom_rec.occupancy,
        atom_rec.temp_fac,
        atom_rec.element,
        atom_rec.charge,
        res)
    # If atom does not exist in the residue, create an Atom
    if !haskey(atoms(res), atom_rec.atom_name)
        res[atom_rec.atom_name] = at
    # Atom exists in the residue, atom names match and alt loc IDs are different
    elseif isa(res[atom_rec.atom_name], Atom) &&
            atom_rec.alt_loc_id != altlocid(res[atom_rec.atom_name])
        # If we are removing disorder and the new atom is preferred to the old one, replace the old one
        if remove_disorder &&
                choosedefaultaltlocid(at, res[atom_rec.atom_name]) == atom_rec.alt_loc_id
            res[atom_rec.atom_name] = at
        # If we are not removing disorder, create a new disordered atom container and add both atoms
        elseif !remove_disorder
            res[atom_rec.atom_name] = DisorderedAtom(Dict(
                atom_rec.alt_loc_id => at,
                altlocid(res[atom_rec.atom_name]) => res[atom_rec.atom_name]
            ), choosedefaultaltlocid(at, res[atom_rec.atom_name]))
        end
    # A disordered atom container already exists and the alt loc ID is not taken
    elseif isa(res[atom_rec.atom_name], DisorderedAtom) &&
            !(atom_rec.alt_loc_id in altlocids(res[atom_rec.atom_name]))
        # Add the new atom to the disordered atom container
        res[atom_rec.atom_name][atom_rec.alt_loc_id] = at
        # If the default alt loc requires changing, change it
        if choosedefaultaltlocid(defaultatom(res[atom_rec.atom_name]), atom) != defaultaltlocid(res[atom_rec.atom_name])
            res[atom_rec.atom_name] = DisorderedAtom(
                        res[atom_rec.atom_name],
                        atom_rec.alt_loc_id)
        end
    else
        error("Two copies of the same atom have the same alternative location ID:\n" *
              "$(res[atom_rec.atom_name])$(atom_rec)")
    end
end


# Internal function to form ordered sub-element lists after parsing
function fixlists!(struc::ProteinStructure)
    for mod in struc
        for ch in mod
            append!(ch.res_list, map(resid, sort(collect(values(residues(ch))))))
            for res in ch
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

function fixlists!(res::Residue)
    append!(res.atom_list, map(fullatomname, sort(collect(values(atoms(res))))))
end

fullatomname(at::AbstractAtom) = atomname(at, spaces=true)


"""
Determine which of two `Atom`s representing a disorered atom better qualifies
as the default location. The `Atom` with the highest occupancy is chosen; in the
case of ties the `Atom` with the lowest alternative location ID in alphabetical
order is chosen.
"""
function choosedefaultaltlocid(at_one::Atom, at_two::Atom)
    if occupancy(at_one) > occupancy(at_two) ||
            (occupancy(at_one) == occupancy(at_two) &&
            Int(altlocid(at_one)) < Int(altlocid(at_two)))
        return altlocid(at_one)
    else
        return altlocid(at_two)
    end
end


"""
Determines if an `AbstractAtom` represents a non-hetero atom, e.g. came from an
ATOM record in a Protein Data Bank (PDB) file.
"""
stdatomselector(at::AbstractAtom) = !ishetatom(at)


"""
Determines if an `AbstractAtom` is a hetero atom, e.g. came from a HETATM
record in a Protein Data Bank (PDB) file.
"""
hetatomselector(at::AbstractAtom) = ishetatom(at)


"""
Determines if an `AbstractAtom` has its atom name in the given `Set` or
`Vector`.
`spaces` determines whether surrounding whitespace is retained
(default `false`).
"""
function atomnameselector(at::AbstractAtom,
                    atom_names::Set{String};
                    spaces::Bool=false)
    return atomname(at, spaces=spaces) in atom_names
end

# Set is faster but Vector method is retained for ease of use
function atomnameselector(at::AbstractAtom,
                    atom_names::Vector{String};
                    spaces::Bool=false)
    return atomname(at, spaces=spaces) in atom_names
end


"`Set` of C-alpha atom names."
const calphaatomnames = Set(["CA"])


"""
Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
C-alpha atom.
"""
function calphaselector(at::AbstractAtom)
    return stdatomselector(at) && atomnameselector(at, calphaatomnames)
end


"`Set` of C-beta atom names."
const cbetaatomnames = Set(["CB"])


"""
Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
C-beta atom, or a C-alpha atom in glycine.
"""
function cbetaselector(at::AbstractAtom)
    return stdatomselector(at) &&
        (atomnameselector(at, cbetaatomnames) ||
        (resname(at) == "GLY" && atomnameselector(at, calphaatomnames)))
end


"`Set` of protein backbone atom names."
const backboneatomnames = Set(["CA", "N", "C"])


"""
Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
protein backbone atom.
"""
function backboneselector(at::AbstractAtom)
    return stdatomselector(at) && atomnameselector(at, backboneatomnames)
end


"""
Determines if an `AbstractAtom` corresponds to a heavy (non-hydrogen) atom and
is not a hetero-atom.
"""
heavyatomselector(at::AbstractAtom) = stdatomselector(at) && !hydrogenselector(at)


"""
Determines if an `AbstractResidue` or `AbstractAtom` has its residue name in
the given `Set` or `Vector`.
"""
function resnameselector(el::Union{AbstractResidue, AbstractAtom},
                    res_names::Set{String})
    return resname(el) in res_names
end

function resnameselector(el::Union{AbstractResidue, AbstractAtom},
                    res_names::Vector{String})
    return resname(el) in res_names
end


"`Set` of residue names corresponding to water."
const waterresnames = Set(["HOH"])


"""
Determines if an `AbstractResidue` or `AbstractAtom` represents a water
molecule.
"""
function waterselector(el::Union{AbstractResidue, AbstractAtom})
    return resnameselector(el, waterresnames)
end


"""
Determines if an `AbstractResidue` represents a standard protein residue,
e.g. consists of ATOM records in a Protein Data Bank (PDB) file.
"""
stdresselector(res::AbstractResidue) = !ishetres(res)


"""
Determines if an `AbstractResidue` represents a hetero molecule, e.g. consists
of HETATM records in a Protein Data Bank (PDB) file.
"""
hetresselector(res::AbstractResidue) = ishetres(res)


"""
Determines whether an `AbstractAtom` or `AbstractResidue` is disordered, i.e.
has multiple locations in the case of atoms or multiple residue names (point
mutants) in the case of residues.
"""
disorderselector(at::AbstractAtom) = isdisorderedatom(at)

disorderselector(res::AbstractResidue) = isdisorderedres(res)


# Either the element is H or the element field is empty, the atom name contains
#   an H and there are no letters in the atom name before H
# For example atom names "1H" and "H1" would be hydrogens but "NH1" would not
"""
Determines if an `AbstractAtom` represents hydrogen. Uses the element field
where possible, otherwise uses the atom name.
"""
function hydrogenselector(at::AbstractAtom)
    return element(at, spaces=false) == "H" ||
        (element(at, spaces=false) == "" &&
        'H' in atomname(at) &&
        !ismatch(r"[a-zA-Z]", atomname(at)[1:findfirst(atomname(at), 'H')-1]))
end


# Residues are ordered in a chain with all hetero residues at the end
# For obtaining sequence we will re-order them numerically
function AminoAcidSequence(ch::Chain, residue_selectors::Function...)
    return AminoAcidSequence(
        sort(collectresidues(ch, residue_selectors...), by=resnumber))
end

function AminoAcidSequence{T <: AbstractResidue}(res::Vector{T})
    seq = AminoAcid[]
    for i in 1:length(res)
        if resname(res[i]) in keys(threeletter_to_aa)
            push!(seq, threeletter_to_aa[resname(res[i])])
        else
            push!(seq, AA_X)
        end
        if i+1 <= length(res) && resnumber(res[i+1]) - resnumber(res[i]) > 1
            append!(seq, [AA_Gap for _ in 1:(resnumber(res[i+1]) - resnumber(res[i]) - 1)])
        end
    end
    return AminoAcidSequence(seq)
end


# Descriptive showing of elements

function Base.show(io::IO, struc::ProteinStructure)
    if countmodels(struc) > 0
        mod = defaultmodel(struc)
        println(io, "Name                        -  ", structurename(struc))
        println(io, "Number of models            -  ", countmodels(struc))
        println(io, "Chain(s)                    -  ", join(chainids(mod)))
        println(io, "Number of residues          -  ", countresidues(mod, stdresselector))
        println(io, "Number of point mutations   -  ", countresidues(mod, stdresselector, disorderselector))
        println(io, "Number of other molecules   -  ", countresidues(mod, hetresselector) - countresidues(mod, hetresselector, waterselector))
        println(io, "Number of water molecules   -  ", countresidues(mod, hetresselector, waterselector))
        println(io, "Number of atoms             -  ", countatoms(mod, stdatomselector))
        println(io, "Number of hydrogens         -  ", countatoms(mod, stdatomselector, hydrogenselector))
        println(io, "Number of disordered atoms  -  ", countatoms(mod, stdatomselector, disorderselector))
    else
        println(io, "Name                        -  ", structurename(struc))
        println(io, "Number of models            -  0")
    end
end

function Base.show(io::IO, mod::Model)
    println(io, "Model number                -  ", modelnumber(mod))
    println(io, "Chain(s)                    -  ", join(chainids(mod)))
    println(io, "Number of residues          -  ", countresidues(mod, stdresselector))
    println(io, "Number of point mutations   -  ", countresidues(mod, stdresselector, disorderselector))
    println(io, "Number of other molecules   -  ", countresidues(mod, hetresselector) - countresidues(mod, hetresselector, waterselector))
    println(io, "Number of water molecules   -  ", countresidues(mod, hetresselector, waterselector))
    println(io, "Number of atoms             -  ", countatoms(mod, stdatomselector))
    println(io, "Number of hydrogens         -  ", countatoms(mod, stdatomselector, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(mod, stdatomselector, disorderselector))
end

function Base.show(io::IO, ch::Chain)
    println(io, "Chain ID                    -  ", chainid(ch))
    println(io, "Number of residues          -  ", countresidues(ch, stdresselector))
    println(io, "Number of point mutations   -  ", countresidues(ch, stdresselector, disorderselector))
    println(io, "Number of other molecules   -  ", countresidues(ch, hetresselector) - countresidues(ch, hetresselector, waterselector))
    println(io, "Number of water molecules   -  ", countresidues(ch, hetresselector, waterselector))
    println(io, "Number of atoms             -  ", countatoms(ch, stdatomselector))
    println(io, "Number of hydrogens         -  ", countatoms(ch, stdatomselector, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(ch, stdatomselector, disorderselector))
end

function Base.show(io::IO, res::Residue)
    println(io, "Residue ID                  -  ", resid(res, full=true))
    println(io, "Residue name                -  ", resname(res))
    println(io, "Number of atoms             -  ", countatoms(res))
    println(io, "Number of hydrogens         -  ", countatoms(res, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(res, disorderselector))
end

function Base.show(io::IO, dis_res::DisorderedResidue)
    println(io, "Residue ID                  -  ", resid(dis_res, full=true))
    for res_name in resnames(dis_res)
        println(io, "Residue name                -  ", res_name)
        println(io, "Number of atoms             -  ", countatoms(disorderedres(dis_res, res_name)))
        println(io, "Number of hydrogens         -  ", countatoms(disorderedres(dis_res, res_name), hydrogenselector))
        println(io, "Number of disordered atoms  -  ", countatoms(disorderedres(dis_res, res_name), disorderselector))
    end
end

Base.show(io::IO, at::Atom) = println(io, pdbline(at)...)

Base.showcompact(io::IO, at::Atom) = print(io, pdbline(at)...)

function Base.show(io::IO, dis_at::DisorderedAtom)
    for at in dis_at
        show(io, at)
    end
end

Base.show(io::IO, at_rec::AtomRecord) = println(io, pdbline(at_rec)...)

Base.showcompact(io::IO, at_rec::AtomRecord) = print(io, pdbline(at_rec)...)
