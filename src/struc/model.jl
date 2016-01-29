export AbstractAtom,
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
    applyselector,
    applyselector!,
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


"""An atom represented in a PDB file - either an `Atom` or a `DisorderedAtom`."""
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

"""A container to hold different versions off the same atom."""
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
    atoms::Dict{AbstractString, AbstractAtom}
end

Residue(name::AbstractString, chain_id::Char, number::Int, ins_code::Char, het_res::Bool) = Residue(name, chain_id, number, ins_code, het_res, Dict())
Residue(name::AbstractString, chain_id::Char, number::Int, ins_code::Char) = Residue(name, chain_id, number, ins_code, false, Dict())
Residue(name::AbstractString, chain_id::Char, number::Int) = Residue(name, chain_id, number, ' ', false, Dict())

"""A container to hold different versions of the same residue (point mutations)"""
immutable DisorderedResidue <: AbstractResidue
    names::Dict{AbstractString, Residue}
    default::AbstractString
end


"""A chain from a PDB file."""
immutable Chain <: StrucElement
    id::Char
    residues::Dict{AbstractString, AbstractResidue}
end

Chain(id::Char) = Chain(id, Dict())


"""A model from a PDB file."""
immutable Model <: StrucElement
    number::Int
    chains::Dict{Char, Chain}
end

Model(number::Int) = Model(number, Dict())


"""A container for multiple models from a PDB file"""
immutable Structure <: StrucElement
    name::AbstractString
    models::Dict{Int, Model}
end

Structure(name::AbstractString) = Structure(name, Dict())


typealias AtomList Array{AbstractAtom,1}
typealias ResidueList Array{AbstractResidue,1}
typealias ChainList Array{Chain,1}
typealias ModelList Array{Model,1}

typealias StrucElementOrList Union{StrucElement, AtomList, ResidueList, ChainList, ModelList}


# Allow accessing sub elements contained in an element like a dictionary
# e.g. allows you to do res[atom_name] rather than having to do res.atoms[atom_name]

getindex(disordered_atom::DisorderedAtom, alt_loc_id::Char) = disordered_atom.alt_loc_ids[alt_loc_id]
function setindex!(disordered_atom::DisorderedAtom, atom::Atom, alt_loc_id::Char)
    disordered_atom.alt_loc_ids[alt_loc_id] = atom
end

getindex(res::Residue, atom_name::AbstractString) = res.atoms[atom_name]
function setindex!(res::Residue, atom::AbstractAtom, atom_name::AbstractString)
    res.atoms[atom_name] = atom
end

getindex(disordered_res::DisorderedResidue, res_name::AbstractString) = disordered_res.names[res_name]
function setindex!(disordered_res::DisorderedResidue, res::Residue, res_name::AbstractString)
    disordered_res.names[res_name] = res
end

getindex(chain::Chain, res_id::AbstractString) = chain.residues[res_id]
function setindex!(chain::Chain, res::AbstractResidue, res_id::AbstractString)
    chain.residues[res_id] = res
end

getindex(chain::Chain, res_no::Int) = chain.residues[string(res_no)]
function setindex!(chain::Chain, res::AbstractResidue, res_no::Int)
    chain.residues[string(res_no)] = res
end

getindex(model::Model, chain_id::Char) = model.chains[chain_id]
function setindex!(model::Model, chain::Chain, chain_id::Char)
    model.chains[chain_id] = chain
end

getindex(struc::Structure, model_no::Int) = struc.models[model_no]
function setindex!(struc::Structure, model::Model, model_no::Int)
    struc.models[model_no] = model
end

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

setx!(atom::Atom, x::Real) = atom.coords[1] = x
sety!(atom::Atom, y::Real) = atom.coords[2] = y
setz!(atom::Atom, z::Real) = atom.coords[3] = z

function setcoords!(atom::Atom, coords::Array{Float64,1})
    @assert length(coords) == 3 "3 coordinates must be given"
    setx!(atom, coords[1])
    sety!(atom, coords[2])
    setz!(atom, coords[3])
end

settempfac!(atom::Atom, temp_fac::Real) = atom.temp_fac = temp_fac

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


function getresid(element::Union{AbstractResidue, AbstractAtom}; full::Bool=false)
    res_id = strip("$(getresnumber(element))$(getinscode(element))")
    ishetero(element) ? res_id = "H_$res_id"
    full ? res_id = "$(res_id):$(getchainid(element))"
    return res_id
end


getresname(res::Residue) = res.name
getchainid(res::Residue) = res.chain_id
getresnumber(res::Residue) = res.number
getinscode(res::Residue) = res.ins_code
ishetres(res::Residue) = res.het_res
getatomnames(res::Residue) = sort(collect(keys(res.atoms)), by= atom_name -> getserial(res[atom_name]))

getdefaultresname(disordered_res::DisorderedResidue) = disordered_res.default
getdefaultresidue(disordered_res::DisorderedResidue) = disordered_res[getdefaultresname(disordered_res)]
getresnames(disordered_res::DisorderedResidue) = sort(collect(keys(disordered_res.res_names))) # Is alphabetical the correct sorting?

getresname(disordered_res::DisorderedResidue) = getdefaultresname(disordered_res)
getchainid(disordered_res::DisorderedResidue) = getchainid(getdefaultresidue(disordered_res))
getresnumber(disordered_res::DisorderedResidue) = getresnumber(getdefaultresidue(disordered_res))
getinscode(disordered_res::DisorderedResidue) = getinscode(getdefaultresidue(disordered_res))
ishetres(disordered_res::DisorderedResidue) = ishetres(getdefaultresidue(disordered_res))
geatomnames(disordered_res::DisorderedResidue) = getatomnames(getdefaultresidue(disordered_res))

ishetero(res::AbstractResidue) = ishetres(res)

function setdefaultresname!(disordered_res::DisorderedResidue, res_name::AbstractString)
    @assert res_name in getresnames(disordered_res) "The new default residue name must be present in the residue"
    disordered_res = DisorderedResidue(disordered_res.res_names, res_name)
end


getchainid(chain::Chain) = chain.id

getresids(chain::Chain) = sortresids(collect(keys(chain.residues)))

function sortresids(res_ids::Array{ASCIIString,1})
    sorted_res_ids = sort(res_ids, by= res_id -> chain[res_id].ins_code)
    sort!(sorted_res_ids, by= res_id -> chain[res_id].res_no)
    sort!(sorted_res_ids, by= res_id -> chain[res_id].het_res)
    return sorted_res_ids
end


getmodelnumber(model::Model) = model.number
getchainids(model::Model) = sort(collect(keys(model.chains)))


getstrucname(struc::Structure) = struc.name
getmodelnumbers(struc::Structure) = sort(collect(keys(struc.models)))

setstrucname!(struc::Structure, name::AbstractString) = struc.name = name


# Iterators to yield sub elements when looping over an element

start(::Structure) = 1
next(struc::Structure, state) = (struc[getmodelnumbers(struc)[state]], state + 1)
done(struc::Structure, state) = state > length(getmodelnumbers(struc))
eltype(::Type{Structure}) = Model
length(struc::Structure) = length(getmodelnumbers(struc))

start(::Model) = 1
next(model::Model, state) = (model[getchainids(model)[state]], state + 1)
done(model::Model, state) = state > length(getchainids(model))
eltype(::Type{Model}) = Chain
length(model::Model) = length(getchainids(model))

start(::Chain) = 1
next(chain::Chain, state) = (chain[getresids(chain)[state]], state + 1)
done(chain::Chain, state) = state > length(getresids(chain))
eltype(::Type{Chain}) = AbstractResidue
length(chain::Chain) = length(getresids(chain))

start(::Residue) = 1
next(res::Residue, state) = (res[getatomnames(res)[state]], state + 1)
done(res::Residue, state) = state > length(getatomnames(res))
eltype(::Type{Residue}) = AbstractAtom
length(res::Residue) = length(getatomnames(res))

start(::DisorderedResidue) = 1
next(disordered_res::DisorderedResidue, state) = (getdefaultresidue(disordered_res)[getatomnames(disordered_res)[state]], state + 1)
done(disordered_res::DisorderedResidue, state) = state > length(getatomnames(disordered_res))
eltype(::Type{DisorderedResidue}) = AbstractAtom
length(disordered_res::DisorderedResidue) = length(getatomnames(disordered_res))

start(::Atom) = 1
next(atom::Atom, state) = (atom, state + 1)
done(atom::Atom, state) = state > 1
eltype(::Type{Atom}) = Atom
length(atom::Atom) = 1

start(::DisorderedAtom) = 1
next(disordered_atom::DisorderedAtom, state) = (disordered_atom[getaltlocids(disordered_atom)[state]], state + 1)
done(disordered_atom::DisorderedAtom, state) = state > length(getaltlocids(disordered_atom))
eltype(::Type{DisorderedAtom}) = Atom
length(disordered_atom::DisorderedAtom) = length(getaltlocids(disordered_atom))


# Is deepcopy needed?
function applyselector(element_list::Union{ResidueList, AtomList}, args...)
    new_list = copy(element_list)
    applyselector!(new_list, args...)
    return new_list
end

function applyselector!(element_list::Union{ResidueList, AtomList}, args...)
    for selector_function in args
        filter!(element_list, element -> selector_function(element))
    end
end


collect(element::Union{Chain, AbstractResidue, AbstractAtom}, args...) = applyselector(collect(element), args...)

# Map here
collectresidues(struc::Structure, args...) = collectresidues(struc[1], args...)
collectresidues(model::Model, args...) = [collectresidues(chain, args...) for chain in model]
collectresidues(chain::Chain, args...) = collect(chain, args...)
collectresidues(res::AbstractResidue, args...) = applyselector([res], args...)
collectresidues(atom::AbstractAtom, args...) = applyselector(organise(atom), args...)
collectresidues(models::ModelList, args...) = [collectresidues(model, args...) for model in models]
collectresidues(chains::ChainList, args...) = [collectresidues(chain, args...) for chain in chains]
collectresidues(residues::ResidueList, args...) = applyselector(residues, args...)
collectresidues(atoms::AtomList, args...) = applyselector(organise(atoms), args...)

collectatoms(struc::Structure, args...) = collectatoms(struc[1], args...)
collectatoms(model::Model, args...) = [collectatoms(chain, args...) for chain in model]
collectatoms(chain::Chain, args...) = [collectatoms(res, args...) for res in chain]
collectatoms(res::AbstractResidue, args...) = collect(res, args...)
collectatoms(atom::AbstractAtom, args...) = collect(atom, args...)
collectatoms(models::ModelList, args...) = [collectatoms(model, args...) for model in models]
collectatoms(chains::ChainList, args...) = [collectatoms(chain, args...) for chain in chains]
collectatoms(residues::ResidueList, args...) = [collectatoms(res, args...) for res in residues]
collectatoms(atoms::AtomList, args...) = [collectatoms(atom, args...) for atom in atoms]


countmodels(struc::Structure) = length(struc)
countchains(struc::Structure) = length(struc[1])
countchains(chain::Chain) = length(chain)
countresidues(element::StrucElementOrList, args...) = length(collectresidues(element, args...))
countatoms(element::StrucElementOrList, args...) = length(collectatoms(element, args...))


# Returns true if the number of DisorderedAtoms or DisorderedResidues is greater than 0
isdisorderd(element::StrucElementOrList) = countatoms(element, disorderselector) > 0 || countresidues(element, disorderselector) > 0


function organise(models::ModelList; struc_name::AbstractString="")
    struc = Structure(struc_name)
    for model in models
        @assert !(getmodelnumber(model) in struc) "Multiple models with the same model number found - cannot organise into a structure"
        struc[getmodelnumber(model)] = model
    end
    return struc
end


function organise(chains::ChainList; model_number::Int=1)
    model = Model(model_number)
    for chain in chains
        @assert !(getchainid(chain) in model) "Multiple chains with the same chain ID found - cannot organise into a model"
        model[getchainid(chain)] = chain
    end
    return model
end


function organise(residues::ResidueList)
    chains = Dict{Char, Dict{AbstractString, AbstractResidue}}()
    for res in residue
        chain_id = getchainid(res)
        if !(chain_id in chains)
            chains[chain_id] = Dict{AbstractString, AbstractResidue}()
        end
        res_id = getresid(res)
        @assert !(res_id in chains[chain_id]) "Multiple residues with the same residue ID found - cannot organise into chains"
        chains[chain_id][res_id] = res
    end
    return [Chain(chain_id, chains[chain_id]) for chain_id in sort(collect(keys(chains)))]
end


function organise(atoms::AtomList)
    residues = Dict{AbstractString, AbstractResidue}()
    for atom in atoms
        res_id = getresid(atom)
        # If the residue doesn't exist, create it
        if !(res_id in residues)
            residues[res_id] = Residue()
            current_residue = residues[res_id]
        # The disordered residue container could already exist
        elseif typeof(residues[res_id]) == DisorderedResidue
            # If the res name isn't present, need to add a new Residue
            if !(getresname(atom) in getresnames(residues[res_id]))
                residues[res_id][getresname(atom)] = Residue(
                    getresname(atom),
                    getchainid(atom),
                    getresnumber(atom),
                    getinscode(atom),
                    ishetatom(atom)
                )
            end
            current_residue = residues[res_id][getresname(atom)]
        # The disorered residue container doesn't exist and needs creating
        elseif getresname(atom) != getresname(residues[res_id])
            # The default res name is the first added
            residues[res_id] = DisorderedResidue(Dict(
                getresname(residues[res_id]) => residues[res_id],
                getresname(atom) => Residue(
                    getresname(atom),
                    getchainid(atom),
                    getresnumber(atom),
                    getinscode(atom),
                    ishetatom(atom)
                )
            ), getresname(residues[res_id]))
            current_residue = residues[res_id][getresname(atom)]
        # Otherwise the residue exists with the same residue name and we can try to add to it
        else
            current_residue = residues[res_id]
        end
        @assert !(getatomname(atom) in getatomnames(current_residue)) "Multiple atoms with the same atom name on the same residue - cannot organise into residues"
        current_residue[getatomname(atom)] = atom
    end
    return [residues[res_id] for res_id in sortresids(collect(keys(residues)))]
end


# This looks slow - instead take a whole list and do it at once like above?
function addtoatomlist!(atoms::AtomList, new_atom::Atom, remove_disorder::Bool=false)
    new_res_id = getresid(new_atom)
    new_atom_name = getatomname(new_atom)
    already_added = false
    # Find if there is an existing atom with the same name on the same residue
    for (i, atom) in enumerate(atoms)
        if getresid(atom) == new_res_id && getatomname(atom) == new_atom_name
            new_alt_loc_id = getaltlocid(new_atom)
            # A disordered atom container already exists and the alt loc ID is not taken
            if typeof(atom) == DisorderedAtom && !(new_alt_loc_id in getaltlocids(atom))
                # Add the new atom to the disordered atom container
                atom[new_alt_loc_id] = new_atom
                # If the default alt loc requires changing, change it
                choosedefaultaltlocid(getdefaultatom(atom), new_atom) != getdefaultaltlocid(atom) ? setdefaultaltlocid!(atom, new_alt_loc_id)
            # The atom already exists and the alt loc ID is not taken
            elseif typeof(atom) == Atom && new_alt_loc_id != getaltlocid(atom)
                # If we are removing disorder and the new atom is preferred to the old one, replace the old one
                if remove_disorder && choosedefaultaltlocid(atom, new_atom) == new_alt_loc_id
                    atom = new_atom
                # If we are not removing disorder, create a new disordered atom container
                else
                    atom = DisorderedAtom(Dict(getaltlocid(atom) => atom, new_alt_loc_id => new_atom), choosedefaultaltlocid(atom, new_atom))
                end
            else
                error("Two copies of the same atom have the same alternative location ID")
            end
            already_added = true
            break
        end
    end
    !already_added ? push!(atoms, new_atom)
end


function choosedefaultaltlocid(atom_one::Atom, atom_two::Atom)
    if getoccupancy(atom_one) > getoccupancy(atom_two) ||
            (getoccupancy(atom_one) == getoccupancy(atom_two) &&
            Int(getaltlocid(atom_one)) < Int(getaltlocid(atom_two)))
        return getaltlocid(atom_one)
    else
        return getaltlocid(atom_two)
    end
end


organise(model::Model) = organise([model])
organise(chain::Chain) = organise([chain])
organise(res::AbstractResidue) = organise([res])
organise(atom::AbstractAtom, args...) = organise([atom], args...)

# These can possibly be put as a Union
organisemodel(chains::ChainList; model_number::Int=1) = organise(chains; model_number=model_number)
organisemodel(residues::ResidueList; model_number::Int=1) = organise(organise(residues); model_number=model_number)
organisemodel(atoms::AtomList; model_number::Int=1) = organise(organise(organise(atoms)); model_number=model_number)
organisemodel(chain::Chain; model_number::Int=1) = organisemodel([chain]; model_number=model_number)
organisemodel(res::AbstractResidue; model_number::Int=1) = organisemodel([res]; model_number=model_number)
organisemodel(atom::AbstractAtom; model_number::Int=1) = organisemodel([atom]; model_number=model_number)

organisestruc(models::ModelList; struc_name::AbstractString="") = organise(models; struc_name=struc_name)
organisestruc(chains::ChainList; struc_name::AbstractString="", model_number::Int=1) = organise(organisemodel(chains; model_number=model_number); struc_name=struc_name)
organisestruc(residues::ResidueList; struc_name::AbstractString="", model_number::Int=1) = organise(organisemodel(residues; model_number=model_number); struc_name=struc_name)
organisestruc(atoms::AtomList; struc_name::AbstractString="", model_number::Int=1) = organise(organisemodel(atoms; model_number=model_number); struc_name=struc_name)
organisestruc(model::Model; struc_name::AbstractString="") = organisestruc([model]; struc_name=struc_name)
organisestruc(chain::Chain; struc_name::AbstractString="", model_number::Int=1) = organisestruc([chain]; struc_name=struc_name, model_number=model_number)
organisestruc(res::AbstractResidue; struc_name::AbstractString="", model_number::Int=1) = organisestruc([res]; struc_name=struc_name, model_number=model_number)
organisestruc(atom::AbstractAtom; struc_name::AbstractString="", model_number::Int=1) = organisestruc([atom]; struc_name=struc_name, model_number=model_number)


stdatomselector(atom::AbstractAtom) = !ishetatom(atom)

hetatomselector(atom::AbstractAtom) = ishetatom(atom)

atomnameselector(atom::AbstractAtom, atom_names::Array{AbstractString,1}) = getatomname(atom) in atom_names

const calpha_atom_names = ["CA"]

calphaselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, calpha_atom_names)

const backbone_atom_names = ["CA", "N", "C"]

backboneselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, backbone_atom_names)

const heavy_atom_names = ["C", "CA", "CB", "CD", "CD1", "CD2", "CE", "CE1",
                    "CE2", "CE3", "CG", "CG1", "CG2", "CH2", "CZ", "CZ2",
                    "CZ3", "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2",
                    "NZ", "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH",
                    "SG"]

heavyatomselector(atom::AbstractAtom) = stdatomselector(atom) && atomnameselector(atom, heavy_atom_names)

resnameselector(element::Union{AbstractResidue, AbstractAtom}, res_names::Array{AbstractString,1}) = getresname(element) in res_names

# Any more?
const water_res_names = ["HOH"]

waterselector(element::Union{AbstractResidue, AbstractAtom}) = resnameselector(element, water_res_names)

stdresselector(res::AbstractResidue) = !ishetres(res)

hetresselector(res::AbstractResidue) = ishetres(res)

disorderselector(atom::AbstractAtom) = typeof(atom) == DisorderedAtom

disorderselector(res::AbstractResidue) = typeof(res) == DisorderedResidue

# Either the element is H or the element field is empty, the atom name contains an H and there are no letters in the atom name before H
# For example atom names 1H and H1 would be hydrogens but NH1 would not
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
