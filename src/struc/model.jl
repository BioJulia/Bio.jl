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
    getaltloc,
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
    getdefaultaltloc,
    getdefaultatom,
    getaltlocs,
    ishetero,
    ishetres,
    getatomnames,
    getresids,
    getmodelnumber,
    getchainids,
    getstrucname,
    getmodelnumbers,
    filteratoms,
    filteratoms!,
    unfold,
    unfoldresidues,
    countresidues,
    unfoldatoms,
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


import Base: getindex, setindex!, start, next, done, show


abstract StrucElement


abstract AbstractAtom <: StrucElement

immutable Atom <: AbstractAtom
    het_atom::Bool
    serial::Int
    name::AbstractString
    alt_loc::Char
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

immutable DisorderedAtom <: AbstractAtom
    alt_locs::Dict{Char, Atom}
    default::Char
end


abstract AbstractResidue <: StrucElement

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

immutable DisorderedResidue <: AbstractResidue
    names::Dict{AbstractString, Residue}
    default::AbstractString
end


immutable Chain <: StrucElement
    id::Char
    residues::Dict{AbstractString, AbstractResidue}
end

Chain(id::Char) = Chain(id, Dict())


immutable Model <: StrucElement
    number::Int
    chains::Dict{Char, Chain}
end

Model(number::Int) = Model(number, Dict())


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


getindex(disordered_atom::DisorderedAtom, alt_loc::Char) = disordered_atom.alt_locs[alt_loc]

function setindex!(disordered_atom::DisorderedAtom, atom::Atom, alt_loc::Char)
    disordered_atom.alt_locs[alt_loc] = atom
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


ishetatom(atom::Atom) = atom.het_atom
getserial(atom::Atom) = atom.serial
getatomname(atom::Atom) = atom.name
getaltloc(atom::Atom) = atom.alt_loc
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

getdefaultaltloc(disordered_atom::DisorderedAtom) = disordered_atom.default

function setdefaultaltloc!(disordered_atom::DisorderedAtom, alt_loc::Char)
    disordered_atom = DisorderedAtom(disordered_atom.alt_locs, alt_loc)
end

getdefaultatom(disordered_atom::DisorderedAtom) = disordered_atom[getdefaultaltloc(disordered_atom)]
getaltlocs(disordered_atom::DisorderedAtom) = sort(collect(keys(disordered_atom.alt_locs)))

ishetatom(disordered_atom::DisorderedAtom) = ishetatom(getdefaultatom(disordered_atom))
getserial(disordered_atom::DisorderedAtom) = getserial(getdefaultatom(disordered_atom))
getatomname(disordered_atom::DisorderedAtom) = getatomname(getdefaultatom(disordered_atom))
getaltloc(disordered_atom::DisorderedAtom) = getdefaultaltloc(disordered_atom)
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


getresname(res::Residue) = res.name
getchainid(res::Residue) = res.chain_id
getresnumber(res::Residue) = res.number
getinscode(res::Residue) = res.ins_code
ishetres(res::Residue) = res.het_res
getatomnames(res::Residue) = sort(collect(keys(res.atoms)), by= atom_name -> getserial(res[atom_name]))


function getresid(element::Union{AbstractResidue, AbstractAtom}; full::Bool=false)
    res_id = strip("$(getresnumber(element))$(getinscode(element))")
    if ishetero(element)
        res_id = "H_$res_id"
    end
    if full
        return "$(res_id):$(getchainid(element))"
    else
        return res_id
    end
end

getdefaultresname(disordered_res::DisorderedResidue) = disordered_res.default
# Is immutatable
#setdefaultresname!(disordered_res::DisorderedResidue, res_name::AbstractString)
getdefaultresidue(disordered_res::DisorderedResidue) = disordered_res[getdefaultresname(disordered_res)]
getresnames(disordered_res::DisorderedResidue) = sort(collect(keys(disordered_res.res_names)))

getresname(disordered_res::DisorderedResidue) = getdefaultresname(disordered_res)
getchainid(disordered_res::DisorderedResidue) = getchainid(getdefaultresidue(disordered_res))
getresnumber(disordered_res::DisorderedResidue) = getresnumber(getdefaultresidue(disordered_res))
getinscode(disordered_res::DisorderedResidue) = getinscode(getdefaultresidue(disordered_res))
ishetres(disordered_res::DisorderedResidue) = ishetres(getdefaultresidue(disordered_res))
geatomnames(disordered_res::DisorderedResidue) = getatomnames(getdefaultresidue(disordered_res))

ishetero(res::AbstractResidue) = ishetres(res)


getchainid(chain::Chain) = chain.id

function getresids(chain::Chain)
    res_ids = sort(collect(keys(chain.residues)), by= res_id -> chain[res_id].ins_code)
    sort!(res_ids, by= res_id -> chain[res_id].res_no)
    sort!(res_ids, by= res_id -> chain[res_id].het_res)
    return res_ids
end


getmodelnumber(model::Model) = model.number
getchainids(model::Model) = sort(collect(keys(model.chains)))


getstrucname(struc::Structure) = struc.name
getmodelnumbers(struc::Structure) = sort(collect(keys(struc.models)))
countmodels(struc::Structure) = length(struc.models)



start(::Structure) = 1
next(struc::Structure, state) = (struc[getmodelnumbers(struc)[state]], state + 1)
done(struc::Structure, state) = state > length(getmodelnumbers(struc))

start(::Model) = 1
next(model::Model, state) = (model[getchainids(model)[state]], state + 1)
done(model::Model, state) = state > length(getchainids(model))

start(::Chain) = 1
next(chain::Chain, state) = (residues[state], state + 1)
done(chain::Chain, state) = state > length(getresids(chain))

start(::Residue) = 1
next(res::Residue, state) = (res[getatomnames(res)[state]], state + 1)
done(res::Residue, state) = state > length(getatomnames(res))

start(::DisorderedResidue) = 1
next(disordered_res::DisorderedResidue, state) = (getdefaultresidue(disordered_res)[getatomnames(disordered_res)[state]], state + 1)
done(disordered_res::DisorderedResidue, state) = state > length(getatomnames(disordered_res))

start(::Atom) = 1
next(atom::Atom, state) = (atom, state + 1)
done(atom::Atom, state) = state > 1

start(::DisorderedAtom) = 1
next(disordered_atom::DisorderedAtom, state) = (disordered_atom[getaltlocs(disordered_atom)[state]], state + 1)
done(disordered_atom::DisorderedAtom, state) = state > length(getaltlocs(disordered_atom))


function filteratoms(atoms::AtomList, args...)
    new_atoms = copy(atoms)
    filteratoms!(new_atoms, args...)
    return new_atoms
end


function filteratoms!(atoms::AtomList, args...)
    for selection_function in args
        filter!(atoms, atom -> selection_function(atom))
    end
end


unfold(element::StrucElement) = [sub_element for sub_element in element]
unfold(atom_container::Union{AbstractResidue, AbstractAtom}, args...) = filteratoms(unfold(atom_container), args...)

unfoldresidues(struc::Structure) = [unfoldresidues(model) for model in struc]
unfoldresidues(model::Model) = [unfoldresidues(chain) for chain in model]
unfoldresidues(chain::Chain) = unfold(chain)
unfoldresidues(res::AbstractResidue) = [res]
unfoldresidues(atom::AbstractAtom) = organise(atom)
unfoldresidues(models::ModelList) = [unfoldresidues(model) for model in models]
unfoldresidues(chains::ChainList) = [unfoldresidues(chain) for chain in chains]
unfoldresidues(residues::ResidueList) = residues
unfoldresidues(atoms::AtomList) = organise(atoms)

countresidues(element::StrucElementOrList) = length(unfoldresidues(element))

unfoldatoms(struc::Structure, args...) = [unfoldatoms(model, args...) for model in struc]
unfoldatoms(model::Model, args...) = [unfoldatoms(chain, args...) for chain in model]
unfoldatoms(chain::Chain, args...) = [unfoldatoms(res, args...) for res in chain]
unfoldatoms(res::AbstractResidue, args...) = unfold(res, args...)
unfoldatoms(atom::AbstractAtom, args...) = unfold(atom, args...)
unfoldatoms(models::ModelList, args...) = [unfoldatoms(model, args...) for model in models]
unfoldatoms(chains::ChainList, args...) = [unfoldatoms(chain, args...) for chain in chains]
unfoldatoms(residues::ResidueList, args...) = [unfoldatoms(res, args...) for res in residues]
unfoldatoms(atoms::AtomList, args...) = [unfoldatoms(atom, args...) for atom in atoms]

countatoms(element::StrucElementOrList, args...) = length(unfoldatoms(element, args...))


function organise(models::ModelList; struc_name::AbstractString="")
    model_dict = Dict{Int, Model}()
    for (model_number, model) in models
        @assert !(model_number in model_dict) "Multiple models with the same model number found - cannot organise into a structure"
        model_dict[model_number] = model
    end
    return Struc(struc_name, model_dict)
end


function organise(chains::ChainList; model_number::Int=1)
    chain_dict = Dict{Char, Chain}()
    for (chain_id, chain) in chains
        @assert !(chain_id in chain_dict) "Multiple chains with the same chain ID found - cannot organise into a model"
        chain_dict[chain_id] = chain
    end
    return Model(model_number, chain_dict)
end


# Currently assumes disordered residues are already sorted (also atoms)
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


# Should this organise into DisorderedAtoms where appropriate?
# This needs to deal with disordered res
function organise(atoms::AtomList)
    residues = Dict{AbstractString, Dict{AbstractString, AbstractAtom}}()
    for atom in atoms
        res_id = getresid(atom)
        if !(res_id in residues)
            residues[res_id] = Dict{AbstractString, AbstractAtom}()
        end
        atom_name = getatomname(atom)
        @assert !(atom_name in residues[res_id]) "Multiple atoms with the same atom name on the same residue - cannot organise into residues"
        residues[res_id][atom_name] = atom
    end
    return
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

calphaselector(atom::AbstractAtom) = atomnameselector(atom, calpha_atom_names)

const backbone_atom_names = ["CA", "N", "C"]

backboneselector(atom::AbstractAtom) = atomnameselector(atom, backbone_atom_names)

const heavy_atom_names = ["C", "CA", "CB", "CD", "CD1", "CD2", "CE", "CE1",
                    "CE2", "CE3", "CG", "CG1", "CG2", "CH2", "CZ", "CZ2",
                    "CZ3", "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2",
                    "NZ", "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH",
                    "SG"]

heavyatomselector(atom::AbstractAtom) = atomnameselector(atom, heavy_atom_names)

resnameselector(atom::AbstractAtom, res_names::Array{AbstractString,1}) = getresname(atom) in res_names

const water_res_names = []

waterselector(atom::AbstractAtom) = resnameselector(atom, water_res_names)


function show(io::IO, struc::Structure)
    model = struc[1]
    println(io, "Name         -  ", getstrucname(struc))
    println(io, "No models    -  ", countmodels(struc))
    println(io, "Chain(s)     -  ", join(getchains(model)))
    println(io, "No residues  -  ", countres(model, het_atom=false))
    println(io, "No atoms     -  ", countatoms(model, het_atom=false))
    println(io, "No hetatoms  -  ", countatoms(model, std_atom=false))
end


function show(io::IO, model::Model)
    println(io, "Model number -  ", model.number)
    println(io, "Chain(s)     -  ", join(getchainids(model)))
    println(io, "No residues  -  ", countresidues(model, stdresselector))
    println(io, "No atoms     -  ", countatoms(model, hetatomselector))
    println(io, "No hetatoms  -  ", countatoms(model, stdatomselector))
end


function show(io::IO, chain::Chain)
    println(io, "Chain ID     -  ", getchainid(chain))
    println(io, "No residues  -  ", countresidues(chain, stdresselector))
    println(io, "No atoms     -  ", countatoms(model, hetatomselector))
    println(io, "No hetatoms  -  ", countatoms(model, stdatomselector))
end


function show(io::IO, res::Residue)
    println(io, "Res name     -  ", getresname(res))
    println(io, "Res ID       -  ", getresid(res, full=true))
    println(io, "No atoms     -  ", countatoms(model, hetatomselector))
    println(io, "No hetatoms  -  ", countatoms(model, stdatomselector))
end


function show(io::IO, atom::Atom)
    print(io, getpdbline(atom)...)
end
