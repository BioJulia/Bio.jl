export StructuralElement,
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
    backbone_atom_names,
    backboneselector,
    heavy_atom_names,
    heavyatomselector,
    resnameselector,
    water_res_names,
    waterselector,
    stdresselector,
    hetresselector,
    disorderselector,
    hydrogenselector


"""A protein structural element."""
abstract StructuralElement


"""An atom represented in a PDB file - either an `Atom` or a
`DisorderedAtom`."""
abstract AbstractAtom <: StructuralElement

"""An atom record from a PDB file."""
immutable Atom <: AbstractAtom
    het_atom::Bool
    serial::Int
    name::ASCIIString
    alt_loc_id::Char
    res_name::ASCIIString
    chain_id::Char
    res_number::Int
    ins_code::Char
    coords::Vector{Float64}
    occupancy::Float64
    temp_fac::Float64
    element::ASCIIString
    charge::ASCIIString
end

"""A container to hold different versions of the same atom."""
immutable DisorderedAtom <: AbstractAtom
    alt_loc_ids::Dict{Char, Atom}
    default::Char
end


"""A residue (amino acid) or other molecule from a PDB file."""
abstract AbstractResidue <: StructuralElement

"""A residue (amino acid) or other molecule from a PDB file."""
immutable Residue <: AbstractResidue
    name::ASCIIString
    chain_id::Char
    number::Int
    ins_code::Char
    het_res::Bool # Does the residue consist of hetatoms?
    atom_list::Vector{ASCIIString}
    atoms::Dict{ASCIIString, AbstractAtom}
end

Residue(name::ASCIIString, chain_id::Char, number::Int, ins_code::Char, het_res::Bool) = Residue(name, chain_id, number, ins_code, het_res, [], Dict())
# Constructor from a list of atoms
Residue{T <: AbstractAtom}(atoms::Vector{T}) = Residue(resname(atoms[1]), chainid(atoms[1]), resnumber(atoms[1]), inscode(atoms[1]), ishetero(atoms[1]), map(atomname, sort(atoms)), [atomname(atom) => atom for atom in atoms])

"""A container to hold different versions of the same residue
(point mutations)."""
immutable DisorderedResidue <: AbstractResidue
    names::Dict{ASCIIString, Residue}
    default::ASCIIString
end


"""A chain from a PDB file."""
immutable Chain <: StructuralElement
    id::Char
    res_list::Vector{ASCIIString}
    residues::Dict{ASCIIString, AbstractResidue}
end

Chain(id::Char) = Chain(id, [], Dict())


"""A model from a PDB file."""
immutable Model <: StructuralElement
    number::Int
    chains::Dict{Char, Chain}
end

Model(number::Int) = Model(number, Dict())
Model() = Model(1)


"""A container for multiple models from a PDB file."""
immutable ProteinStructure <: StructuralElement
    name::ASCIIString
    models::Dict{Int, Model}
end

ProteinStructure(name::ASCIIString) = ProteinStructure(name, Dict())
ProteinStructure() = ProteinStructure("")


"""A protein structural element or list of structural elements."""
typealias StructuralElementOrList Union{StructuralElement, Vector{AbstractAtom}, Vector{Atom}, Vector{DisorderedAtom}, Vector{AbstractResidue}, Vector{Residue}, Vector{DisorderedResidue}, Vector{Chain}, Vector{Model}}


# Allow accessing sub elements contained in an element like a dictionary
# e.g. allows you to do res[atom_name] rather than res.atoms[atom_name]
# setindex! should be used with caution as it is not checked and can lead to inconsistencies

# Accessing a DisorderedAtom with a character returns the Atom with that alt loc ID
Base.getindex(disordered_atom::DisorderedAtom, alt_loc_id::Char) = disordered_atom.alt_loc_ids[alt_loc_id]
function Base.setindex!(disordered_atom::DisorderedAtom, atom::Atom, alt_loc_id::Char)
    disordered_atom.alt_loc_ids[alt_loc_id] = atom
end

# Accessing a Residue with an ASCIIString returns the AbstractAtom with that atom name
Base.getindex(res::Residue, atom_name::ASCIIString) = res.atoms[atom_name]
function Base.setindex!(res::Residue, atom::AbstractAtom, atom_name::ASCIIString)
    res.atoms[atom_name] = atom
end

# Accessing a DisorderedResidue with an ASCIIString returns the AbstractAtom in the default Residue with that atom name
# This is not necessarily intuitive, it may be expected to return the Residue with that residue name
# However this way accessing an AbstractResidue always returns an AbstractAtom
Base.getindex(disordered_res::DisorderedResidue, atom_name::ASCIIString) = disordered_res.names[defaultresname(disordered_res)][atom_name]
function Base.setindex!(disordered_res::DisorderedResidue, atom::AbstractAtom, atom_name::ASCIIString)
    disordered_res.names[defaultresname(disordered_res)][atom_name] = atom
end

# Accessing a Chain with an ASCIIString returns the AbstractResidue with that residue ID
Base.getindex(chain::Chain, res_id::ASCIIString) = chain.residues[res_id]
function Base.setindex!(chain::Chain, res::AbstractResidue, res_id::ASCIIString)
    chain.residues[res_id] = res
end

# Accessing a Chain with an Int returns the AbstractResidue with that residue ID converted to a string
Base.getindex(chain::Chain, res_number::Int) = chain.residues[string(res_number)]
function Base.setindex!(chain::Chain, res::AbstractResidue, res_number::Int)
    chain.residues[string(res_number)] = res
end

# Accessing a Model with a Char returns the Chain with that chain ID
Base.getindex(model::Model, chain_id::Char) = model.chains[chain_id]
function Base.setindex!(model::Model, chain::Chain, chain_id::Char)
    model.chains[chain_id] = chain
end

# Accessing a ProteinStructure with an Int returns the Model with that model number
Base.getindex(struc::ProteinStructure, model_number::Int) = struc.models[model_number]
function Base.setindex!(struc::ProteinStructure, model::Model, model_number::Int)
    struc.models[model_number] = model
end

# Accessing a ProteinStructure with a Char returns the Chain with that chain ID on the default model
Base.getindex(struc::ProteinStructure, chain_id::Char) = defaultmodel(struc)[chain_id]
function Base.setindex!(struc::ProteinStructure, chain::Chain, chain_id::Char)
    defaultmodel(struc)[chain_id] = chain
end


# Getters and setters for structural elements

# Atom getters/setters
ishetatom(atom::Atom) = atom.het_atom
serial(atom::Atom) = atom.serial
atomname(atom::Atom) = atom.name
altlocid(atom::Atom) = atom.alt_loc_id
resname(atom::Atom) = atom.res_name
chainid(atom::Atom) = atom.chain_id
resnumber(atom::Atom) = atom.res_number
inscode(atom::Atom) = atom.ins_code
x(atom::Atom) = atom.coords[1]
y(atom::Atom) = atom.coords[2]
z(atom::Atom) = atom.coords[3]
coords(atom::Atom) = atom.coords
occupancy(atom::Atom) = atom.occupancy
tempfac(atom::Atom) = atom.temp_fac
element(atom::Atom) = atom.element
charge(atom::Atom) = atom.charge

ishetero(atom::AbstractAtom) = ishetatom(atom)
isdisorderedatom(::Atom) = false
atomid(atom::Atom) = (resid(atom; full=true), resname(atom), atomname(atom))

x!(atom::Atom, x::Real) = (atom.coords[1] = x; nothing)
y!(atom::Atom, y::Real) = (atom.coords[2] = y; nothing)
z!(atom::Atom, z::Real) = (atom.coords[3] = z; nothing)

function coords!(atom::Atom, coords::Vector{Float64})
    @assert length(coords) == 3 "3 coordinates must be given"
    x!(atom, coords[1])
    y!(atom, coords[2])
    z!(atom, coords[3])
    return nothing
end

# DisorderedAtom getters/setters
defaultaltlocid(disordered_atom::DisorderedAtom) = disordered_atom.default
defaultatom(disordered_atom::DisorderedAtom) = disordered_atom[defaultaltlocid(disordered_atom)]
altlocids(disordered_atom::DisorderedAtom) = sort(collect(keys(disordered_atom.alt_loc_ids)), by= alt_loc_id -> serial(disordered_atom[alt_loc_id]))

ishetatom(disordered_atom::DisorderedAtom) = ishetatom(defaultatom(disordered_atom))
serial(disordered_atom::DisorderedAtom) = serial(defaultatom(disordered_atom))
atomname(disordered_atom::DisorderedAtom) = atomname(defaultatom(disordered_atom))
altlocid(disordered_atom::DisorderedAtom) = defaultaltlocid(disordered_atom)
resname(disordered_atom::DisorderedAtom) = resname(defaultatom(disordered_atom))
chainid(disordered_atom::DisorderedAtom) = chainid(defaultatom(disordered_atom))
resnumber(disordered_atom::DisorderedAtom) = resnumber(defaultatom(disordered_atom))
inscode(disordered_atom::DisorderedAtom) = inscode(defaultatom(disordered_atom))
x(disordered_atom::DisorderedAtom) = x(defaultatom(disordered_atom))
y(disordered_atom::DisorderedAtom) = y(defaultatom(disordered_atom))
z(disordered_atom::DisorderedAtom) = z(defaultatom(disordered_atom))
# Returns a Vector of atomic coordinates - compare with coordarray which returns a 2D array
coords(disordered_atom::DisorderedAtom) = coords(defaultatom(disordered_atom))
occupancy(disordered_atom::DisorderedAtom) = occupancy(defaultatom(disordered_atom))
tempfac(disordered_atom::DisorderedAtom) = tempfac(defaultatom(disordered_atom))
element(disordered_atom::DisorderedAtom) = element(defaultatom(disordered_atom))
charge(disordered_atom::DisorderedAtom) = charge(defaultatom(disordered_atom))

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

function resid(element::Union{AbstractResidue, AbstractAtom}; full::Bool=false)
    if ishetero(element)
        if full
            inscode(element) == ' ' ? "H_$(resnumber(element)):$(chainid(element))" : "H_$(resnumber(element))$(inscode(element)):$(chainid(element))"
        else
            inscode(element) == ' ' ? "H_$(resnumber(element))" : "H_$(resnumber(element))$(inscode(element))"
        end
    else
        if full
            inscode(element) == ' ' ? "$(resnumber(element)):$(chainid(element))" : "$(resnumber(element))$(inscode(element)):$(chainid(element))"
        else
            inscode(element) == ' ' ? "$(resnumber(element))" : "$(resnumber(element))$(inscode(element))"
        end
    end
end


# Residue getters/setters
resname(res::Residue) = res.name
chainid(res::Residue) = res.chain_id
resnumber(res::Residue) = res.number
inscode(res::Residue) = res.ins_code
ishetres(res::Residue) = res.het_res
atomnames(res::Residue) = res.atom_list
atoms(res::Residue) = res.atoms

ishetero(res::AbstractResidue) = ishetres(res)
isdisorderedres(::Residue) = false

# DisorderedResidue getters/setters
disorderedres(disordered_res::DisorderedResidue, res_name::ASCIIString) = disordered_res.names[res_name]
defaultresname(disordered_res::DisorderedResidue) = disordered_res.default
defaultresidue(disordered_res::DisorderedResidue) = disordered_res.names[defaultresname(disordered_res)]
# Orders default residue name first then others alphabetically
resnames(disordered_res::DisorderedResidue) = sort(collect(keys(disordered_res.names)), lt= (res_name_one, res_name_two) -> (isless(res_name_one, res_name_two) && res_name_two != defaultresname(disordered_res)) || res_name_one == defaultresname(disordered_res))

resname(disordered_res::DisorderedResidue) = defaultresname(disordered_res)
chainid(disordered_res::DisorderedResidue) = chainid(defaultresidue(disordered_res))
resnumber(disordered_res::DisorderedResidue) = resnumber(defaultresidue(disordered_res))
inscode(disordered_res::DisorderedResidue) = inscode(defaultresidue(disordered_res))
ishetres(disordered_res::DisorderedResidue) = ishetres(defaultresidue(disordered_res))
atomnames(disordered_res::DisorderedResidue) = atomnames(defaultresidue(disordered_res))
atoms(disordered_res::DisorderedResidue) = atoms(defaultresidue(disordered_res))

isdisorderedres(::DisorderedResidue) = true

# Constructor acts as a setter for the default residue name
function DisorderedResidue(disordered_res::DisorderedResidue, default::ASCIIString)
    @assert default in resnames(disordered_res) "The new default residue name must be present in the residue"
    return DisorderedResidue(disordered_res.names, default)
end


# Chain getters/setters
chainid(chain::Chain) = chain.id
resids(chain::Chain) = chain.res_list
residues(chain::Chain) = chain.residues


# Model getters/setters
modelnumber(model::Model) = model.number
chainids(model::Model) = sortchainids(collect(keys(model.chains)))
chains(model::Model) = model.chains


# ProteinStructure getters/setters
structurename(struc::ProteinStructure) = struc.name
modelnumbers(struc::ProteinStructure) = sort(collect(keys(struc.models)))
models(struc::ProteinStructure) = struc.models

defaultmodel(struc::ProteinStructure) = struc.models[modelnumbers(struc)[1]]
chainids(struc::ProteinStructure) = countmodels(struc) > 0 ? chainids(defaultmodel(struc)) : Char[]


# Explicit sorting functions

Base.sort!{T <: AbstractAtom}(atoms::Vector{T}) = sort!(atoms, by=serial)

function Base.sort{T <: AbstractAtom}(atoms::Vector{T})
    atoms_copy = copy(atoms)
    sort!(atoms_copy)
    return atoms_copy
end


# Assumes all on same chain
function Base.sort!{T <: AbstractResidue}(residues::Vector{T})
    sort!(residues, by=inscode)
    sort!(residues, by=resnumber)
    sort!(residues, by=ishetres)
    sort!(residues, by=chainid, lt=chainidisless)
end

function Base.sort{T <: AbstractResidue}(residues::Vector{T})
    residues_copy = copy(residues)
    sort!(residues_copy)
    return residues_copy
end


# Chains are ordered by character sorting except the empty chain ID comes last
sortchainids(chain_ids::Vector{Char}) = sort(chain_ids, lt=chainidisless)

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


"""Returns a copy of a `Vector{AbstractResidue}` or `Vector{AbstractAtom}` with all elements that do not
satisfy `selector_functions...` removed."""
function applyselectors{T <: Union{AbstractResidue, AbstractAtom}}(element_list::Vector{T}, selector_functions::Function...)
    new_list = copy(element_list)
    applyselectors!(new_list, selector_functions...)
    return new_list
end

"""Runs `applyselectors` in place."""
function applyselectors!{T <: Union{AbstractResidue, AbstractAtom}}(element_list::Vector{T}, selector_functions::Function...)
    for selector_function in selector_functions
        filter!(selector_function, element_list)
    end
end

applyselectors{T <: Union{AbstractResidue, AbstractAtom}}(element_list::Vector{T}) = copy(element_list)
applyselectors!{T <: Union{AbstractResidue, AbstractAtom}}(element_list::Vector{T}) = element_list


"""Returns a `Vector{AbstractResidue}` of the residues in an element.
Only elements that satisfy `selector_functions...` are returned."""
collectresidues(struc::ProteinStructure, selector_functions::Function...) = countmodels(struc) > 0 ? collectresidues(defaultmodel(struc), selector_functions...) : AbstractResidue[]
collectresidues(chain::Chain, selector_functions::Function...) = applyselectors(collect(chain), selector_functions...)
collectresidues(res::AbstractResidue, selector_functions::Function...) = applyselectors(AbstractResidue[res], selector_functions...)
collectresidues(atom::AbstractAtom, selector_functions::Function...) = applyselectors(organise(atom), selector_functions...)
# Note output is always Vector{AbstractResidue} unless input was Vector{Residue} or Vector{DisorderedResidue}, in which case output is same type as input type
collectresidues{T <: AbstractResidue}(residues::Vector{T}, selector_functions::Function...) = sort(applyselectors(residues, selector_functions...))
collectresidues{T <: AbstractAtom}(atoms::Vector{T}, selector_functions::Function...) = applyselectors(organise(atoms), selector_functions...)

function collectresidues(element::Union{Model, Vector{Model}, Vector{Chain}}, selector_functions::Function...)
    residues = AbstractResidue[]
    for sub_element in element
        append!(residues, collectresidues(sub_element, selector_functions...))
    end
    return residues
end


"""Return an `Vector{AbstractAtom}` of the atoms in an element, ordered by serial.
Only elements that satisfy `selector_functions...` are returned."""
collectatoms(struc::ProteinStructure, selector_functions::Function...) = countmodels(struc) > 0 ? collectatoms(defaultmodel(struc), selector_functions...) : AbstractAtom[]
collectatoms(res::AbstractResidue, selector_functions::Function...) = applyselectors(collect(res), selector_functions...)
collectatoms(atom::AbstractAtom, selector_functions::Function...) = applyselectors(AbstractAtom[atom], selector_functions...)
# Note output is always Vector{AbstractAtom} unless input was Vector{Atom} or Vector{DisorderedAtom}, in which case output is same type as input type
collectatoms{T <: AbstractAtom}(atoms::Vector{T}, selector_functions::Function...) = sort(applyselectors(atoms, selector_functions...))

function collectatoms(element::Union{Model, Chain, Vector{Model}, Vector{Chain}, Vector{AbstractResidue}, Vector{Residue}, Vector{DisorderedResidue}}, selector_functions::Function...)
    atoms = AbstractAtom[]
    for sub_element in element
        append!(atoms, collectatoms(sub_element, selector_functions...))
    end
    return atoms
end


countmodels(struc::ProteinStructure) = length(struc)
countchains(struc::ProteinStructure) = countmodels(struc) > 0 ? length(defaultmodel(struc)) : 0
countchains(model::Model) = length(model)
countresidues(element::StructuralElementOrList, selector_functions::Function...) = length(collectresidues(element, selector_functions...))
countatoms(element::StructuralElementOrList, selector_functions::Function...) = length(collectatoms(element, selector_functions...))


"""Organise a `StructuralElementOrList` into the next level up the heirarchy. An
`Vector{AbstractAtom}` becomes a `Vector{AbstractResidue}`, a `Vector{AbstractResidue}` becomes a `Vector{Chain}`, a
`Vector{Chain}` becomes a `Model` and a `Vector{Model}` becomes a `ProteinStructure`."""
function organise(models::Vector{Model}; structure_name::ASCIIString="")
    # Organise a Vector{Model} into a ProteinStructure
    struc = ProteinStructure(structure_name)
    for model in models
        @assert !(modelnumber(model) in modelnumbers(struc)) "Multiple models with the same model number found - cannot organise into a protein structure"
        struc[modelnumber(model)] = model
    end
    return struc
end


# Organise a Vector{Chain} into a Model
function organise(chains::Vector{Chain}; model_number::Int=1)
    model = Model(model_number)
    for chain in chains
        @assert !(chainid(chain) in chainids(model)) "Multiple chains with the same chain ID found - cannot organise into a model"
        model[chainid(chain)] = chain
    end
    return model
end


# Organise a Vector{AbstractResidue} into a Vector{Chain}
function organise{T <: AbstractResidue}(residues::Vector{T})
    chains = Dict{Char, Dict{ASCIIString, AbstractResidue}}()
    for res in residues
        chain_id = chainid(res)
        !haskey(chains, chain_id) ? chains[chain_id] = Dict{ASCIIString, AbstractResidue}() : nothing
        res_id = resid(res)
        @assert !haskey(chains[chain_id], res_id) "Multiple residues with the same residue ID found - cannot organise into chains"
        chains[chain_id][res_id] = res
    end
    return [Chain(chain_id, map(resid, sort(collect(values(chains[chain_id])))), chains[chain_id]) for chain_id in sortchainids(collect(keys(chains)))]
end


# Organise a Vector{AbstractAtom} into a Vector{AbstractResidue}
function organise{T <: AbstractAtom}(atoms::Vector{T})
    # Key is chain ID, value is Dict where key is residue ID and value is list of atoms
    residues = Dict{Char, Dict{ASCIIString, Vector{AbstractAtom}}}()
    # Key is chain ID, value is Dict where key is residue ID, value is Dict where key is residue name and value is list of atoms
    disordered_residues = Dict{Char, Dict{ASCIIString, Dict{ASCIIString, Vector{AbstractAtom}}}}()
    # Key is chain ID, value is Dict where key is residue ID and value is default residue name
    defaults = Dict{Char, Dict{ASCIIString, ASCIIString}}()
    for atom in atoms
        # Create chain Dict if required
        chain_id = chainid(atom)
        if !haskey(residues, chain_id)
            residues[chain_id] = Dict{ASCIIString, Vector{AbstractAtom}}()
            disordered_residues[chain_id] = Dict{ASCIIString, Dict{ASCIIString, Vector{AbstractAtom}}}()
            defaults[chain_id] = Dict{ASCIIString, ASCIIString}()
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
        # Fix quote and atom name check
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


organise(model::Model; structure_name::ASCIIString="") = organise([model]; structure_name=structure_name)
organise(chain::Chain; model_number::Int=1) = organise([chain]; model_number=model_number)
organise(res::AbstractResidue) = organise(AbstractResidue[res])
organise(atom::AbstractAtom) = organise(AbstractAtom[atom])


"""Form a list of `AbstractAtom`s from a list of `Atom`s. Combines disordered
atoms into disordered atom containers unless `remove_disorder` is `true`, in
which case removes all but one location for disordered atoms."""
function formatomlist(atoms::Vector{Atom}; remove_disorder::Bool=false)
    # Key is (residue ID, residue name, atom name)
    atom_dic = Dict{Tuple{ASCIIString, ASCIIString, ASCIIString}, AbstractAtom}()
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
    return sort(collect(values(atom_dic)), by=serial) # Is by serial redundant?
end


"""Determine which of two `Atom`s representing a disorered atom better qualifies
as the default location.
The `Atom` with the highest occupancy is chosen; in the case of ties the `Atom`
with the lowest alternative location ID in alphabetical order is chosen."""
function choosedefaultaltlocid(atom_one::Atom, atom_two::Atom)
    if occupancy(atom_one) > occupancy(atom_two) ||
            (occupancy(atom_one) == occupancy(atom_two) &&
            Int(altlocid(atom_one)) < Int(altlocid(atom_two)))
        return altlocid(atom_one)
    else
        return altlocid(atom_two)
    end
end


# These can possibly be put as a Union
"""Organise elements into a `Model`."""
organisemodel(chains::Vector{Chain}; model_number::Int=1) = organise(chains; model_number=model_number)
organisemodel{T <: AbstractResidue}(residues::Vector{T}; model_number::Int=1) = organise(organise(residues); model_number=model_number)
organisemodel{T <: AbstractAtom}(atoms::Vector{T}; model_number::Int=1) = organise(organise(organise(atoms)); model_number=model_number)
organisemodel(chain::Chain; model_number::Int=1) = organisemodel([chain]; model_number=model_number)
organisemodel(res::AbstractResidue; model_number::Int=1) = organisemodel(AbstractResidue[res]; model_number=model_number)
organisemodel(atom::AbstractAtom; model_number::Int=1) = organisemodel(AbstractAtom[atom]; model_number=model_number)


"""Organise elements into a `ProteinStructure`."""
organisestructure(models::Vector{Model}; structure_name::ASCIIString="") = organise(models; structure_name=structure_name)
organisestructure(chains::Vector{Chain}; structure_name::ASCIIString="", model_number::Int=1) = organise(organisemodel(chains; model_number=model_number); structure_name=structure_name)
organisestructure{T <: AbstractResidue}(residues::Vector{T}; structure_name::ASCIIString="", model_number::Int=1) = organise(organisemodel(residues; model_number=model_number); structure_name=structure_name)
organisestructure{T <: AbstractAtom}(atoms::Vector{T}; structure_name::ASCIIString="", model_number::Int=1) = organise(organisemodel(atoms; model_number=model_number); structure_name=structure_name)
organisestructure(model::Model; structure_name::ASCIIString="") = organisestructure([model]; structure_name=structure_name)
organisestructure(chain::Chain; structure_name::ASCIIString="", model_number::Int=1) = organisestructure([chain]; structure_name=structure_name, model_number=model_number)
organisestructure(res::AbstractResidue; structure_name::ASCIIString="", model_number::Int=1) = organisestructure(AbstractResidue[res]; structure_name=structure_name, model_number=model_number)
organisestructure(atom::AbstractAtom; structure_name::ASCIIString="", model_number::Int=1) = organisestructure(AbstractAtom[atom]; structure_name=structure_name, model_number=model_number)


"""Determines if an `AbstractAtom` is a non-hetero atom, i.e. came from an ATOM
record."""
stdatomselector(atom::AbstractAtom) = !ishetatom(atom)

"""Determines if an `AbstractAtom` is a hetero atom, i.e. came from a HETATM
record."""
hetatomselector(atom::AbstractAtom) = ishetatom(atom)

"""Determines if an `AbstractAtom` has its atom name in the given `Array`."""
atomnameselector(atom::AbstractAtom, atom_names::Set{ASCIIString}) = atomname(atom) in atom_names
# Set is faster but Vector method is retained for ease of use
atomnameselector(atom::AbstractAtom, atom_names::Vector{ASCIIString}) = atomname(atom) in atom_names

"""`Array` of C-alpha atom names."""
const calpha_atom_names = Set(["CA"])

"""Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
C-alpha atom."""
calphaselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, calpha_atom_names)

"""`Array` of protein backbone atom names."""
const backbone_atom_names = Set(["CA", "N", "C"])

"""Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
protein backbone atom."""
backboneselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, backbone_atom_names)

"""`Array` of protein standard heavy atom names."""
const heavy_atom_names = Set(["C", "CA", "CB", "CD", "CD1", "CD2",
                    "CE", "CE1", "CE2", "CE3", "CG", "CG1", "CG2", "CH2", "CZ",
                    "CZ2", "CZ3", "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1",
                    "NH2", "NZ", "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1",
                    "OH", "SG"])

"""Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
standard protein heavy atom."""
heavyatomselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, heavy_atom_names)

"""Determines if an `AbstractResidue` or `AbstractAtom` has its resiudue name in
the given `Array`."""
resnameselector(element::Union{AbstractResidue, AbstractAtom}, res_names::Set{ASCIIString}) = resname(element) in res_names
resnameselector(element::Union{AbstractResidue, AbstractAtom}, res_names::Vector{ASCIIString}) = resname(element) in res_names

"""`Array` of residue names corresponding to water."""
const water_res_names = Set(["HOH"])

"""Determines if an `AbstractResidue` or `AbstractAtom` represents a water
molecule."""
waterselector(element::Union{AbstractResidue, AbstractAtom}) = resnameselector(element, water_res_names)

"""Determines if an `AbstractResidue` represents a standard protein residue,
i.e. consists of ATOM records."""
stdresselector(res::AbstractResidue) = !ishetres(res)

"""Determines if an `AbstractResidue` represents a hetero molecule, i.e.
consists of HETATM records."""
hetresselector(res::AbstractResidue) = ishetres(res)

"""Determines whether an `AbstractAtom` or `AbstractResidue` is disordered, i.e.
has multiple locations in the case of atoms or multiple residue names (point
mutants) in the case of residues."""
disorderselector(atom::AbstractAtom) = isdisorderedatom(atom)
disorderselector(res::AbstractResidue) = isdisorderedres(res)

# Either the element is H or the element field is empty, the atom name contains an H and there are no letters in the atom name before H
# For example atom names 1H and H1 would be hydrogens but NH1 would not
"""Determines if an `AbstractAtom` represents hydrogen. Uses the element field
where possible, otherwise uses the atom name."""
hydrogenselector(atom::AbstractAtom) = element(atom) == "H" || (element(atom) == "" && 'H' in atomname(atom) && !ismatch(r"[a-zA-Z]", atomname(atom)[1:findfirst(atomname(atom), 'H')-1]))


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
        println(io, "Number of hydrogens           -  ", countatoms(model, stdatomselector, hydrogenselector))
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
    println(io, "Number of hydrogens           -  ", countatoms(model, stdatomselector, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(model, stdatomselector, disorderselector))
end

function Base.show(io::IO, chain::Chain)
    println(io, "Chain ID                    -  ", chainid(chain))
    println(io, "Number of residues          -  ", countresidues(chain, stdresselector))
    println(io, "Number of point mutations   -  ", countresidues(chain, stdresselector, disorderselector))
    println(io, "Number of other molecules   -  ", countresidues(chain, hetresselector) - countresidues(chain, hetresselector, waterselector))
    println(io, "Number of water molecules   -  ", countresidues(chain, hetresselector, waterselector))
    println(io, "Number of atoms             -  ", countatoms(chain, stdatomselector))
    println(io, "Number of hydrogens           -  ", countatoms(chain, stdatomselector, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(chain, stdatomselector, disorderselector))
end

function Base.show(io::IO, res::Residue)
    println(io, "Residue ID                  -  ", resid(res; full=true))
    println(io, "Residue name                -  ", resname(res))
    println(io, "Number of atoms             -  ", countatoms(res))
    println(io, "Number of hydrogens           -  ", countatoms(res, hydrogenselector))
    println(io, "Number of disordered atoms  -  ", countatoms(res, disorderselector))
end

function Base.show(io::IO, disordered_res::DisorderedResidue)
    println(io, "Residue ID                  -  ", resid(disordered_res; full=true))
    for res_name in resnames(disordered_res)
        println(io, "Residue name                -  ", res_name)
        println(io, "Number of atoms             -  ", countatoms(disorderedres(disordered_res, res_name)))
        println(io, "Number of hydrogens           -  ", countatoms(disorderedres(disordered_res, res_name), hydrogenselector))
        println(io, "Number of disordered atoms  -  ", countatoms(disorderedres(disordered_res, res_name), disorderselector))
    end
end

function Base.show(io::IO, atom::Atom)
    print(io, pdbline(atom)...)
end

function Base.show(io::IO, disordered_atom::DisorderedAtom)
    for atom in disordered_atom
        show(io, atom)
        println()
    end
end
