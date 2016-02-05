export StrucElement,
    AbstractAtom,
    Atom,
    DisorderedAtom,
    AbstractResidue,
    Residue,
    DisorderedResidue,
    Chain,
    Model,
    Structure,
    AtomList,
    ResidueList,
    ChainList,
    ModelList,
    StrucElementOrList,
    ishetatom,
    getserial,
    getatomname,
    getaltlocid,
    getresname,
    getchainid,
    getresnumber,
    getinscode,
    getx,
    gety,
    getz,
    getcoords,
    getoccupancy,
    gettempfac,
    getelement,
    getcharge,
    setx!,
    sety!,
    setz!,
    setcoords!,
    getdefaultaltlocid,
    getdefaultatom,
    getaltlocids,
    ishetero,
    ishetres,
    getatomnames,
    getresid,
    getdefaultresname,
    getdefaultresidue,
    getresids,
    getmodelnumber,
    getchainids,
    getstrucname,
    getmodelnumbers,
    applyselectors,
    applyselectors!,
    collectresidues,
    countresidues,
    collectatoms,
    countatoms,
    organise,
    organisemodel,
    organisestruc,
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
    waterselector


import Base: getindex,
    setindex!,
    start,
    next,
    done,
    eltype,
    length,
    collect,
    show


"""A protein structural element."""
abstract StrucElement


"""An atom represented in a PDB file - either an `Atom` or a
`DisorderedAtom`."""
abstract AbstractAtom <: StrucElement

"""An atom record from a PDB file."""
immutable Atom <: AbstractAtom
    het_atom::Bool
    serial::Int
    name::AbstractString
    alt_loc_id::Char
    res_name::AbstractString
    chain_id::Char
    res_no::Int
    ins_code::Char
    coords::Array{Float64,1}
    occupancy::Float64
    temp_fac::Float64
    element::AbstractString
    charge::AbstractString
end

"""A container to hold different versions of the same atom."""
immutable DisorderedAtom <: AbstractAtom
    alt_loc_ids::Dict{Char, Atom}
    default::Char
end


"""A residue (amino acid) or other molecule from a PDB file."""
abstract AbstractResidue <: StrucElement

"""A residue (amino acid) or other molecule from a PDB file."""
immutable Residue <: AbstractResidue
    name::AbstractString
    chain_id::Char
    number::Int
    ins_code::Char
    het_res::Bool # Does the residue consist of hetatoms?
    atom_list::Array{AbstractString,1}
    atoms::Dict{AbstractString, AbstractAtom}
end

Residue(name::AbstractString, chain_id::Char, number::Int, ins_code::Char, het_res::Bool) = Residue(name, chain_id, number, ins_code, het_res, AbstractString[], Dict())
Residue(atom::AbstractAtom) = Residue(getresname(atom), getchainid(atom), getresnumber(atom), getinscode(atom), ishetatom(atom))

"""A container to hold different versions of the same residue
(point mutations)."""
immutable DisorderedResidue <: AbstractResidue
    names::Dict{AbstractString, Residue}
    default::AbstractString
end


"""A chain from a PDB file."""
immutable Chain <: StrucElement
    id::Char
    res_list::Array{AbstractString,1}
    residues::Dict{AbstractString, AbstractResidue}
end

Chain(id::Char) = Chain(id, Dict())


"""A model from a PDB file."""
immutable Model <: StrucElement
    number::Int
    chains::Dict{Char, Chain}
end

Model(number::Int) = Model(number, Dict())


"""A container for multiple models from a PDB file."""
immutable Structure <: StrucElement
    name::AbstractString
    models::Dict{Int, Model}
end

Structure(name::AbstractString) = Structure(name, Dict())


"""An `Array{AbstractAtom,1}`."""
typealias AtomList Array{AbstractAtom,1}

"""An `Array{AbstractResidue,1}`."""
typealias ResidueList Array{AbstractResidue,1}

"""An `Array{Chain,1}`."""
typealias ChainList Array{Chain,1}

"""An `Array{Model,1}`."""
typealias ModelList Array{Model,1}

"""A protein structural element or list of structural elements."""
typealias StrucElementOrList Union{StrucElement, AtomList, ResidueList, ChainList, ModelList}


# Allow accessing sub elements contained in an element like a dictionary
# e.g. allows you to do res[atom_name] rather than res.atoms[atom_name]

# Accessing a DisorderedAtom with a character returns the Atom with that alt loc ID
getindex(disordered_atom::DisorderedAtom, alt_loc_id::Char) = disordered_atom.alt_loc_ids[alt_loc_id]
function setindex!(disordered_atom::DisorderedAtom, atom::Atom, alt_loc_id::Char)
    disordered_atom.alt_loc_ids[alt_loc_id] = atom
end

# Accessing a Residue with an AbstractString returns the AbstractAtom with that atom name
getindex(res::Residue, atom_name::AbstractString) = res.atoms[atom_name]
function setindex!(res::Residue, atom::AbstractAtom, atom_name::AbstractString)
    res.atoms[atom_name] = atom
end

# Accessing a DisorderedResidue with an AbstractString returns the AbstractAtom in the default Residue with that atom name
# This is not necessarily intuitive, it may be expected to return the Residue with that residue name
# However this way accessing an AbstractResidue always returns an AbstractAtom
getindex(disordered_res::DisorderedResidue, atom_name::AbstractString) = disordered_res.names[disordered_res.default][atom_name]
function setindex!(disordered_res::DisorderedResidue, atom::AbstractAtom, atom_name::AbstractString)
    disordered_res.names[disordered_res.default][atom_name] = atom
end

# Accessing a Chain with an AbstractString returns the AbstractResidue with that residue ID
getindex(chain::Chain, res_id::AbstractString) = chain.residues[res_id]
function setindex!(chain::Chain, res::AbstractResidue, res_id::AbstractString)
    chain.residues[res_id] = res
end

# Accessing a Chain with an Int returns the AbstractResidue with that residue ID converted to a string
getindex(chain::Chain, res_no::Int) = chain.residues[string(res_no)]
function setindex!(chain::Chain, res::AbstractResidue, res_no::Int)
    chain.residues[string(res_no)] = res
end

# Accessing a Model with a Char returns the Chain with that chain ID
getindex(model::Model, chain_id::Char) = model.chains[chain_id]
function setindex!(model::Model, chain::Chain, chain_id::Char)
    model.chains[chain_id] = chain
end

# Accessing a Structure with an Int returns the Model with that model number
getindex(struc::Structure, model_no::Int) = struc.models[model_no]
function setindex!(struc::Structure, model::Model, model_no::Int)
    struc.models[model_no] = model
end

# Accessing a Structure with a Char returns the Chain with that chain ID on model number 1
getindex(struc::Structure, chain_id::Char) = struc.models[1][chain_id]
function setindex!(struc::Structure, chain::Chain, chain_id::Char)
    struc.models[1][chain_id] = chain
end


# Getters and setters for structural elements

ishetatom(atom::Atom) = atom.het_atom
getserial(atom::Atom) = atom.serial
getatomname(atom::Atom) = atom.name
getaltlocid(atom::Atom) = atom.alt_loc_id
getresname(atom::Atom) = atom.res_name
getchainid(atom::Atom) = atom.chain_id
getresnumber(atom::Atom) = atom.res_no
getinscode(atom::Atom) = atom.ins_code
getx(atom::Atom) = atom.coords[1]
gety(atom::Atom) = atom.coords[2]
getz(atom::Atom) = atom.coords[3]
getcoords(atom::Atom) = atom.coords
getoccupancy(atom::Atom) = atom.occupancy
gettempfac(atom::Atom) = atom.temp_fac
getelement(atom::Atom) = atom.element
getcharge(atom::Atom) = atom.charge

setx!(atom::Atom, x::Real) = (atom.coords[1] = x; nothing)
sety!(atom::Atom, y::Real) = (atom.coords[2] = y; nothing)
setz!(atom::Atom, z::Real) = (atom.coords[3] = z; nothing)

# Return nothing here?
function setcoords!(atom::Atom, coords::Array{Float64,1})
    @assert length(coords) == 3 "3 coordinates must be given"
    setx!(atom, coords[1])
    sety!(atom, coords[2])
    setz!(atom, coords[3])
    return nothing
end

getdefaultaltlocid(disordered_atom::DisorderedAtom) = disordered_atom.default
getdefaultatom(disordered_atom::DisorderedAtom) = disordered_atom[getdefaultaltlocid(disordered_atom)]
getaltlocids(disordered_atom::DisorderedAtom) = sort(collect(keys(disordered_atom.alt_loc_ids)), by= alt_loc_id -> getserial(disordered_atom[alt_loc_id]))

ishetatom(disordered_atom::DisorderedAtom) = ishetatom(getdefaultatom(disordered_atom))
getserial(disordered_atom::DisorderedAtom) = getserial(getdefaultatom(disordered_atom))
getatomname(disordered_atom::DisorderedAtom) = getatomname(getdefaultatom(disordered_atom))
getaltlocid(disordered_atom::DisorderedAtom) = getdefaultaltlocid(disordered_atom)
getresname(disordered_atom::DisorderedAtom) = getresname(getdefaultatom(disordered_atom))
getchainid(disordered_atom::DisorderedAtom) = getchainid(getdefaultatom(disordered_atom))
getresnumber(disordered_atom::DisorderedAtom) = getresnumber(getdefaultatom(disordered_atom))
getinscode(disordered_atom::DisorderedAtom) = getinscode(getdefaultatom(disordered_atom))
getx(disordered_atom::DisorderedAtom) = getx(getdefaultatom(disordered_atom))
gety(disordered_atom::DisorderedAtom) = gety(getdefaultatom(disordered_atom))
getz(disordered_atom::DisorderedAtom) = getz(getdefaultatom(disordered_atom))
getcoords(disordered_atom::DisorderedAtom) = getcoords(getdefaultatom(disordered_atom))
getoccupancy(disordered_atom::DisorderedAtom) = getoccupancy(getdefaultatom(disordered_atom))
gettempfac(disordered_atom::DisorderedAtom) = gettempfac(getdefaultatom(disordered_atom))
getelement(disordered_atom::DisorderedAtom) = getelement(getdefaultatom(disordered_atom))
getcharge(disordered_atom::DisorderedAtom) = getcharge(getdefaultatom(disordered_atom))

ishetero(atom::AbstractAtom) = ishetatom(atom)

function setdefaultaltlocid!(disordered_atom::DisorderedAtom, alt_loc_id::Char)
    @assert alt_loc_id in getaltlocids(disordered_atom) "The new default alternative location ID must be present in the atom"
    disordered_atom = DisorderedAtom(disordered_atom.alt_loc_ids, alt_loc_id)
end

# These coordinate setters only set the default atom coordinates
setx!(disordered_atom::DisorderedAtom, x::Real) = setx!(getdefaultatom(disordered_atom), x)
sety!(disordered_atom::DisorderedAtom, y::Real) = sety!(getdefaultatom(disordered_atom), y)
setz!(disordered_atom::DisorderedAtom, z::Real) = setz!(getdefaultatom(disordered_atom), z)
setcoords!(disordered_atom::DisorderedAtom, coords::Array{Float64,1}) = setcoords!(getdefaultatom(disordered_atom), coords)


function getresid(element::Union{AbstractResidue, AbstractAtom}; full::Bool=false)
    res_id = strip("$(getresnumber(element))$(getinscode(element))")
    ishetero(element) ? res_id = "H_$res_id" : nothing
    full ? res_id = "$(res_id):$(getchainid(element))" : nothing
    return res_id
end


getresname(res::Residue) = res.name
getchainid(res::Residue) = res.chain_id
getresnumber(res::Residue) = res.number
getinscode(res::Residue) = res.ins_code
ishetres(res::Residue) = res.het_res
getatomnames(res::Residue) = res.atom_list

getdisorderedres(disordered_res::DisorderedResidue, res_name::AbstractString) = disordered_res.names[res_name]
getdefaultresname(disordered_res::DisorderedResidue) = disordered_res.default
getdefaultresidue(disordered_res::DisorderedResidue) = disordered_res.names[getdefaultresname(disordered_res)]
getresnames(disordered_res::DisorderedResidue) = sort(collect(keys(disordered_res.names)))

getresname(disordered_res::DisorderedResidue) = getdefaultresname(disordered_res)
getchainid(disordered_res::DisorderedResidue) = getchainid(getdefaultresidue(disordered_res))
getresnumber(disordered_res::DisorderedResidue) = getresnumber(getdefaultresidue(disordered_res))
getinscode(disordered_res::DisorderedResidue) = getinscode(getdefaultresidue(disordered_res))
ishetres(disordered_res::DisorderedResidue) = ishetres(getdefaultresidue(disordered_res))
getatomnames(disordered_res::DisorderedResidue) = getatomnames(getdefaultresidue(disordered_res))

ishetero(res::AbstractResidue) = ishetres(res)

function setdefaultresname!(disordered_res::DisorderedResidue, res_name::AbstractString)
    @assert res_name in getresnames(disordered_res) "The new default residue name must be present in the residue"
    disordered_res = DisorderedResidue(disordered_res.names, res_name)
end


getchainid(chain::Chain) = chain.id

getresids(chain::Chain) = chain.res_list

function sortresids(residues::Dict{AbstractString, AbstractResidue})
    sorted_res_ids = sort(collect(keys(residues)), by= res_id -> getinscode(residues[res_id]))
    sort!(sorted_res_ids, by= res_id -> getresnumber(residues[res_id]))
    sort!(sorted_res_ids, by= res_id -> ishetres(residues[res_id]))
    return sorted_res_ids
end


getmodelnumber(model::Model) = model.number
getchainids(model::Model) = sort(collect(keys(model.chains)))


getstrucname(struc::Structure) = struc.name
getmodelnumbers(struc::Structure) = sort(collect(keys(struc.models)))

setstrucname!(struc::Structure, name::AbstractString) = struc.name = name


# Iterators to yield sub elements when looping over an element

# Iterating over a Structure yields Models
length(struc::Structure) = length(getmodelnumbers(struc))
start(::Structure) = 1
next(struc::Structure, state) = (struc[getmodelnumbers(struc)[state]], state + 1)
done(struc::Structure, state) = state > length(struc)
eltype(::Type{Structure}) = Model

# Iterating over a Model yields Chains
length(model::Model) = length(getchainids(model))
start(::Model) = 1
next(model::Model, state) = (model[getchainids(model)[state]], state + 1)
done(model::Model, state) = state > length(model)
eltype(::Type{Model}) = Chain

# Iterating over a Chain yields AbstractResidues
length(chain::Chain) = length(getresids(chain))
start(::Chain) = 1
next(chain::Chain, state) = (chain[getresids(chain)[state]], state + 1)
done(chain::Chain, state) = state > length(chain)
eltype(::Type{Chain}) = AbstractResidue

# Iterating over a Residue yields AbstractAtoms
length(res::Residue) = length(getatomnames(res))
start(::Residue) = 1
next(res::Residue, state) = (res[getatomnames(res)[state]], state + 1)
done(res::Residue, state) = state > length(res)
eltype(::Type{Residue}) = AbstractAtom

# Iterating over a DisorderedResidue yields AbstractAtoms
# This is not necessarily intuitive, it may be expected to yield Residues
# However this way iterating over an AbstractResidue always yields AbstractAtoms
length(disordered_res::DisorderedResidue) = length(getatomnames(disordered_res))
start(::DisorderedResidue) = 1
next(disordered_res::DisorderedResidue, state) = (getdefaultresidue(disordered_res)[getatomnames(disordered_res)[state]], state + 1)
done(disordered_res::DisorderedResidue, state) = state > length(disordered_res)
eltype(::Type{DisorderedResidue}) = AbstractAtom

# Iterating over an Atom returns itself
# This is not necessarily intuitive, it may be expected to not be an iterator
# However this way iterating over an AbstractAtom always yields Atoms
length(atom::Atom) = 1
start(::Atom) = 1
next(atom::Atom, state) = (atom, state + 1)
done(atom::Atom, state) = state > 1
eltype(::Type{Atom}) = Atom

# Iterating over a DisorderedAtom yields Atoms
length(disordered_atom::DisorderedAtom) = length(getaltlocids(disordered_atom))
start(::DisorderedAtom) = 1
next(disordered_atom::DisorderedAtom, state) = (disordered_atom[getaltlocids(disordered_atom)[state]], state + 1)
done(disordered_atom::DisorderedAtom, state) = state > length(disordered_atom)
eltype(::Type{DisorderedAtom}) = Atom


"""Returns a copy of a `ResidueList` or `AtomList` with all elements that do not
satisfy `args...` removed."""
function applyselectors(element_list::Union{ResidueList, AtomList}, args...)
    new_list = copy(element_list)
    applyselectors!(new_list, args...)
    return new_list
end

"""Runs `applyselectors` in place."""
function applyselectors!(element_list::Union{ResidueList, AtomList}, args...)
    for selector_function in args
        filter!(selector_function, element_list)
    end
end

applyselectors(element_list::Union{ResidueList, AtomList}) = element_list
applyselectors!(element_list::Union{ResidueList, AtomList}) = element_list


"""Returns a `ResidueList` of the residues in an element.
Only elements that satisfy `args...` are returned."""
collectresidues(struc::Structure, args...) = collectresidues(struc[1], args...)
collectresidues(chain::Chain, args...) = applyselectors(collect(chain), args...)
collectresidues(res::AbstractResidue, args...) = applyselectors(AbstractResidue[res], args...)
collectresidues(atom::AbstractAtom, args...) = applyselectors(organise(atom), args...)
collectresidues(residues::ResidueList, args...) = applyselectors(residues, args...)
collectresidues(atoms::AtomList, args...) = applyselectors(organise(atoms), args...)

function collectresidues(element::Union{Model, ModelList, ChainList}, args...)
    residues = AbstractResidue[]
    for sub_element in element
        append!(residues, collectresidues(sub_element, args...))
    end
    return residues
end


"""Return an `AtomList` of the atoms in an element.
Only elements that satisfy `args...` are returned."""
collectatoms(struc::Structure, args...) = collectatoms(struc[1], args...)
collectatoms(res::AbstractResidue, args...) = applyselectors(collect(res), args...)
collectatoms(atom::AbstractAtom, args...) = applyselectors(AbstractAtom[atom], args...)
collectatoms(atoms::AtomList, args...) = applyselectors(atoms, args...)

function collectatoms(element::Union{Model, Chain, ModelList, ChainList, ResidueList}, args...)
    atoms = AbstractAtom[]
    for sub_element in element
        append!(atoms, collectatoms(sub_element, args...))
    end
    return atoms
end


countmodels(struc::Structure) = length(struc)
countchains(struc::Structure) = length(struc[1])
countchains(chain::Chain) = length(chain)
countresidues(element::StrucElementOrList, args...) = length(collectresidues(element, args...))
countatoms(element::StrucElementOrList, args...) = length(collectatoms(element, args...))


# Returns true if the number of DisorderedAtoms or DisorderedResidues is greater than 0
# This needs some thinking about as does not necessarily correspond to a non-blank alt loc ID
isdisordered(element::StrucElementOrList) = countatoms(element, disorderselector) > 0 || countresidues(element, disorderselector) > 0


"""Organise a `StrucElementOrList` into the next level up the heirarchy. An
`AtomList` becomes a `ResidueList`, a `ResidueList` becomes a `ChainList`, a
`ChainList` becomes a `Model` and a `ModelList` becomes a `Structure`."""
function organise(models::ModelList; struc_name::AbstractString="")
    # Organise a ModelList into a Structure
    struc = Structure(struc_name)
    for model in models
        @assert !(getmodelnumber(model) in getmodelnumbers(struc)) "Multiple models with the same model number found - cannot organise into a structure"
        struc[getmodelnumber(model)] = model
    end
    return struc
end


# Organise a ChainList into a Model
function organise(chains::ChainList; model_number::Int=1)
    model = Model(model_number)
    for chain in chains
        @assert !(getchainid(chain) in getchainids(model)) "Multiple chains with the same chain ID found - cannot organise into a model"
        model[getchainid(chain)] = chain
    end
    return model
end


# Organise a ResidueList into a ChainList
function organise(residues::ResidueList)
    chains = Dict{Char, Dict{AbstractString, AbstractResidue}}()
    for res in residues
        chain_id = getchainid(res)
        !(chain_id in keys(chains)) ? chains[chain_id] = Dict{AbstractString, AbstractResidue}() : nothing
        res_id = getresid(res)
        @assert !(res_id in keys(chains[chain_id])) "Multiple residues with the same residue ID found - cannot organise into chains"
        chains[chain_id][res_id] = res
    end
    return [Chain(chain_id, sortresids(chains[chain_id]), chains[chain_id]) for chain_id in sort(collect(keys(chains)))]
end


# Organise an AtomList into a ResidueList
# Perhaps store chain as well so can be properly sorted
function organise(atoms::AtomList)
    # Key is residue ID, value is list of atoms
    residues = Dict{AbstractString, AtomList}()
    # Key is residue ID, value is Dict where key is residue name and value is list of atoms
    disordered_residues = Dict{AbstractString, Dict{AbstractString, AtomList}}()
    # Key is residue ID, value is default residue name
    defaults = Dict{AbstractString, AbstractString}()
    for atom in atoms
        res_id = getresid(atom; full=true)
        if !(res_id in keys(residues)) && !(res_id in keys(disordered_residues))
            residues[res_id] = AbstractAtom[]
            current_residue = residues[res_id]
        # The disordered residue container could already exist
        elseif res_id in keys(disordered_residues)
            # If the res name isn't present, need to add a new Residue
            !(getresname(atom) in keys(disordered_residues[res_id])) ? disordered_residues[res_id][getresname(atom)] = AbstractAtom[] : nothing
            current_residue = disordered_residues[res_id][getresname(atom)]
        # The disorered residue container doesn't exist and needs creating
        elseif getresname(atom) != getresname(residues[res_id][1])
            disordered_residues[res_id] = Dict(
                getresname(residues[res_id][1]) => residues[res_id],
                getresname(atom) => AbstractAtom[]
            )
            # The default res name is the first added
            defaults[res_id] = getresname(residues[res_id][1])
            # Remove Residue now we have created a DisorderedResidue
            delete!(residues, res_id)
            current_residue = disordered_residues[res_id][getresname(atom)]
        else
            current_residue = residues[res_id]
        end
        atom_name = getatomname(atom)
        # Fix quote and atom name check
        index_found = 0
        for i in eachindex(current_residue)
            atom_name == getatomname(current_residue[i]) ? (index_found = i; break) : nothing
        end
        @assert index_found == 0 "Multiple atoms with the same atom name on the same residue - cannot organise into residues:\n$(current_residue[index_found])\n$(atom)"
        push!(current_residue, atom)
    end
    residues_out = AbstractResidue[]
    for res_id in keys(residues)
        atoms = sort(residues[res_id], by= getserial)
        push!(residues_out, Residue(
            getresname(atoms[1]),
            getchainid(atoms[1]),
            getresnumber(atoms[1]),
            getinscode(atoms[1]),
            ishetero(atoms[1]),
            map(getatomname, atoms),
            [getatomname(atom) => atom for atom in atoms]
        ))
    end
    for res_id in keys(disordered_residues)
        new_disordered_res = DisorderedResidue(Dict(), defaults[res_id])
        for res_name in keys(disordered_residues[res_id])
            atoms = sort(disordered_residues[res_id][res_name], by= getserial)
            # Hacky setter
            new_disordered_res.names[res_name] = Residue(
                getresname(atoms[1]),
                getchainid(atoms[1]),
                getresnumber(atoms[1]),
                getinscode(atoms[1]),
                ishetero(atoms[1]),
                map(getatomname, atoms),
                [getatomname(atom) => atom for atom in atoms]
            )
        end
        push!(residues_out, new_disordered_res)
    end
    return residues_out
end


organise(model::Model) = organise([model])
organise(chain::Chain) = organise([chain])
organise(res::AbstractResidue) = organise(AbstractResidue[res])
organise(atom::AbstractAtom) = organise(AbstractAtom[atom])


# Should check none of the atoms are DisorderedAtoms already, or enforce Array{Atom,1}
"""Form a list of `AbstractAtom`s from a list of `Atom`s. Combines disordered
atoms into disordered atom containers unless `remove_disorder` is `true`, in
which case removes all but one location for disordered atoms."""
function formatomlist(atoms::AtomList; remove_disorder::Bool=false)
    # Key is (residue ID, residue name, atom name)
    # Remove this to a separate function getatomid?
    atom_dic = Dict{Tuple{AbstractString, AbstractString, AbstractString}, AbstractAtom}()
    for atom in atoms
        #@assert typeof(atom) != DisorderedAtom "There is already a DisorderedAtom container in the AtomList"
        atom_id = (getresid(atom; full=true), getresname(atom), getatomname(atom))
        # The atom does not exist so we can create it
        if !(atom_id in keys(atom_dic))
            atom_dic[atom_id] = atom
        # A disordered atom container already exists and the alt loc ID is not taken
        elseif typeof(atom_dic[atom_id]) == DisorderedAtom && !(getaltlocid(atom) in getaltlocids(atom_dic[atom_id]))
            # Add the new atom to the disordered atom container
            atom_dic[atom_id][getaltlocid(atom)] = atom
            # If the default alt loc requires changing, change it
            choosedefaultaltlocid(getdefaultatom(atom_dic[atom_id]), atom) != getdefaultaltlocid(atom) ? setdefaultaltlocid!(atom_dic[atom_id], getaltlocid(atom)) : nothing
        # The atom already exists and the alt loc IDs are different
        elseif typeof(atom_dic[atom_id]) == Atom && getaltlocid(atom) != getaltlocid(atom_dic[atom_id])
            # If we are removing disorder and the new atom is preferred to the old one, replace the old one
            if remove_disorder && choosedefaultaltlocid(atom, atom_dic[atom_id]) == getaltlocid(atom)
                atom_dic[atom_id] = atom
            # If we are not removing disorder, create a new disordered atom container and add both atoms
            else
                atom_dic[atom_id] = DisorderedAtom(Dict(
                    getaltlocid(atom) => atom,
                    getaltlocid(atom_dic[atom_id]) => atom_dic[atom_id]
                ), choosedefaultaltlocid(atom, atom_dic[atom_id]))
            end
        else
            error("Two copies of the same atom have the same alternative location ID (note names are stripped of whitespace so identical names with different spacing, e.g. ' CA ' and 'CA  ', are not currently supported):\n$(atom_dic[atom_id])\n$(atom)")
        end
    end
    return sort(collect(values(atom_dic)), by=getserial)
end


"""Determine which of two `Atom`s representing a disorered atom better qualifies
as the default location.
The `Atom` with the highest occupancy is chosen; in the case of ties the `Atom`
with the lowest alternative location ID in alphabetical order is chosen."""
function choosedefaultaltlocid(atom_one::Atom, atom_two::Atom)
    if getoccupancy(atom_one) > getoccupancy(atom_two) ||
            (getoccupancy(atom_one) == getoccupancy(atom_two) &&
            Int(getaltlocid(atom_one)) < Int(getaltlocid(atom_two)))
        return getaltlocid(atom_one)
    else
        return getaltlocid(atom_two)
    end
end


# These can possibly be put as a Union
"""Organise elements into a `Model`."""
organisemodel(chains::ChainList; model_number::Int=1) = organise(chains; model_number=model_number)
organisemodel(residues::ResidueList; model_number::Int=1) = organise(organise(residues); model_number=model_number)
organisemodel(atoms::AtomList; model_number::Int=1) = organise(organise(organise(atoms)); model_number=model_number)
organisemodel(chain::Chain; model_number::Int=1) = organisemodel([chain]; model_number=model_number)
organisemodel(res::AbstractResidue; model_number::Int=1) = organisemodel(AbstractResidue[res]; model_number=model_number)
organisemodel(atom::AbstractAtom; model_number::Int=1) = organisemodel(AbstractAtom[atom]; model_number=model_number)


"""Organise elements into a `Structure`."""
organisestruc(models::ModelList; struc_name::AbstractString="") = organise(models; struc_name=struc_name)
organisestruc(chains::ChainList; struc_name::AbstractString="", model_number::Int=1) = organise(organisemodel(chains; model_number=model_number); struc_name=struc_name)
organisestruc(residues::ResidueList; struc_name::AbstractString="", model_number::Int=1) = organise(organisemodel(residues; model_number=model_number); struc_name=struc_name)
organisestruc(atoms::AtomList; struc_name::AbstractString="", model_number::Int=1) = organise(organisemodel(atoms; model_number=model_number); struc_name=struc_name)
organisestruc(model::Model; struc_name::AbstractString="") = organisestruc([model]; struc_name=struc_name)
organisestruc(chain::Chain; struc_name::AbstractString="", model_number::Int=1) = organisestruc([chain]; struc_name=struc_name, model_number=model_number)
organisestruc(res::AbstractResidue; struc_name::AbstractString="", model_number::Int=1) = organisestruc(AbstractResidue[res]; struc_name=struc_name, model_number=model_number)
organisestruc(atom::AbstractAtom; struc_name::AbstractString="", model_number::Int=1) = organisestruc(AbstractAtom[atom]; struc_name=struc_name, model_number=model_number)


"""Determines if an `AbstractAtom` is a non-hetero atom, i.e. came from an ATOM
record."""
stdatomselector(atom::AbstractAtom) = !ishetatom(atom)

"""Determines if an `AbstractAtom` is a hetero atom, i.e. came from a HETATM
record."""
hetatomselector(atom::AbstractAtom) = ishetatom(atom)

"""Determines if an `AbstractAtom` has its atom name in the given `Array`."""
atomnameselector(atom::AbstractAtom, atom_names::Array{AbstractString,1}) = getatomname(atom) in atom_names

"""`Array` of C-alpha atom names."""
const calpha_atom_names = AbstractString["CA"]

"""Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
C-alpha atom."""
calphaselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, calpha_atom_names)

"""`Array` of protein backbone atom names."""
const backbone_atom_names = AbstractString["CA", "N", "C"]

"""Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
protein backbone atom."""
backboneselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, backbone_atom_names)

"""`Array` of protein standard heavy atom names."""
const heavy_atom_names = AbstractString["C", "CA", "CB", "CD", "CD1", "CD2",
                    "CE", "CE1", "CE2", "CE3", "CG", "CG1", "CG2", "CH2", "CZ",
                    "CZ2", "CZ3", "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1",
                    "NH2", "NZ", "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1",
                    "OH", "SG"]

"""Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
standard protein heavy atom."""
heavyatomselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, heavy_atom_names)

"""Determines if an `AbstractResidue` or `AbstractAtom` has its resiudue name in
the given `Array`."""
resnameselector(element::Union{AbstractResidue, AbstractAtom}, res_names::Array{AbstractString,1}) = getresname(element) in res_names

# Any more?
"""`Array` of residue names corresponding to water."""
const water_res_names = AbstractString["HOH"]

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
disorderselector(atom::AbstractAtom) = typeof(atom) == DisorderedAtom
disorderselector(res::AbstractResidue) = typeof(res) == DisorderedResidue

# Either the element is H or the element field is empty, the atom name contains an H and there are no letters in the atom name before H
# For example atom names 1H and H1 would be hydrogens but NH1 would not
"""Determines if an `AbstractAtom` represents hydrogen. Uses the element field
where possible, otherwise uses the atom name."""
hydrogenselector(atom::AbstractAtom) = getelement(atom) == "H" || (getelement(atom) == "" && 'H' in getatomname(atom) && !ismatch(r"[a-zA-Z]", getatomname(atom)[:findfirst(getatomname(atom), 'H')-1]))


function show(io::IO, struc::Structure)
    model = struc[1]
    println(io, "Name                 -  ", getstrucname(struc))
    println(io, "No models            -  ", countmodels(struc))
    println(io, "Chain(s)             -  ", join(getchainids(model)))
    println(io, "No residues          -  ", countresidues(model, stdresselector))
    println(io, "No point mutations   -  ", countresidues(model, stdresselector, disorderselector))
    println(io, "No other molecules   -  ", countresidues(model, hetresselector) - countresidues(model, hetresselector, waterselector))
    println(io, "No water molecules   -  ", countresidues(model, hetresselector, waterselector))
    println(io, "No atoms             -  ", countatoms(model, stdatomselector))
    println(io, "No H atoms           -  ", countatoms(model, stdatomselector, hydrogenselector))
    println(io, "No disordered atoms  -  ", countatoms(model, stdatomselector, disorderselector))
end

function show(io::IO, model::Model)
    println(io, "Model number         -  ", getmodelnumber(model))
    println(io, "Chain(s)             -  ", join(getchainids(model)))
    println(io, "No residues          -  ", countresidues(model, stdresselector))
    println(io, "No point mutations   -  ", countresidues(model, stdresselector, disorderselector))
    println(io, "No other molecules   -  ", countresidues(model, hetresselector) - countresidues(model, hetresselector, waterselector))
    println(io, "No water molecules   -  ", countresidues(model, hetresselector, waterselector))
    println(io, "No atoms             -  ", countatoms(model, stdatomselector))
    println(io, "No H atoms           -  ", countatoms(model, stdatomselector, hydrogenselector))
    println(io, "No disordered atoms  -  ", countatoms(model, stdatomselector, disorderselector))
end

function show(io::IO, chain::Chain)
    println(io, "Chain ID             -  ", getchainid(chain))
    println(io, "No residues          -  ", countresidues(chain, stdresselector))
    println(io, "No point mutations   -  ", countresidues(chain, stdresselector, disorderselector))
    println(io, "No other molecules   -  ", countresidues(chain, hetresselector) - countresidues(chain, hetresselector, waterselector))
    println(io, "No water molecules   -  ", countresidues(chain, hetresselector, waterselector))
    println(io, "No atoms             -  ", countatoms(chain, stdatomselector))
    println(io, "No H atoms           -  ", countatoms(chain, stdatomselector, hydrogenselector))
    println(io, "No disordered atoms  -  ", countatoms(chain, stdatomselector, disorderselector))
end

function show(io::IO, res::Residue)
    println(io, "Res ID               -  ", getresid(res, full=true))
    println(io, "Res name             -  ", getresname(res))
    println(io, "No atoms             -  ", countatoms(res))
    println(io, "No H atoms           -  ", countatoms(res, hydrogenselector))
    println(io, "No disordered atoms  -  ", countatoms(res, disorderselector))
end

function show(io::IO, disordered_res::DisorderedResidue)
    println(io, "Res ID               -  ", getresid(disordered_res, full=true))
    for res_name in getresnames(disordered_res)
        println(io, "Res name             -  ", res_name)
        println(io, "No atoms             -  ", countatoms(disordered_res[res_name]))
        println(io, "No H atoms           -  ", countatoms(disordered_res[res_name], hydrogenselector))
        println(io, "No disordered atoms  -  ", countatoms(disordered_res[res_name], disorderselector))
    end
end

function show(io::IO, atom::Atom)
    print(io, getpdbline(atom)...)
end

function show(io::IO, disordered_atom::DisorderedAtom)
    for atom in disordered_atom
        show(io, atom)
        println()
    end
end
