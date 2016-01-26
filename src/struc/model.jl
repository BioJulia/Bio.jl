import Base: start, next, done, getindex, setindex!


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
    x::Float64
    y::Float64
    z::Float64
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
    het_res::Bool # Does the residue consist only of hetatoms?
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
getx(atom::Atom) = atom.x
gety(atom::Atom) = atom.y
getz(atom::Atom) = atom.z
getcoords(atom::Atom) = [atom.x, atom.y, atom.z]
getoccupancy(atom::Atom) = atom.occupancy
gettempfac(atom::Atom) = atom.temp_fac
getelement(atom::Atom) = atom.element
getcharge(atom::Atom) = atom.charge

getdefaultaltloc(disordered_atom::DisorderedAtom) = disordered_atom.default

# No, it's immutable...
function setdefaultaltloc!(disordered_atom::DisorderedAtom, alt_loc::Char)
    disordered_atom.default = alt_loc
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


getresname(res::Residue) = res.name
getchainid(res::Residue) = res.chain_id
getresnumber(res::Residue) = res.number
getinscode(res::Residue) = res.ins_code
ishetres(res::Residue) = res.het_res
getatomnames(res::Residue) = sort(collect(keys(res.atoms)), by= atom_name -> getserial(res[atom_name]))

# NB if this is repeated for atoms the hetero tag is ambiguous
function getresid(res::Residue; full::Bool=false)
    res_id = strip("$(getresnumber(res))$(getinscode(res))")
    if ishetres(res)
        res_id = "H_$res_id"
    end
    if full
        return "$(res_id):$(getchainid(res))"
    else
        return res_id
    end
end

getdefaultresname(disordered_res::DisorderedResidue) = disordered_res.default
# Is immutatable
setdefaultresname!(disordered_res::DisorderedResidue, res_name::AbstractString)
getdefaultresidue(disordered_res::DisorderedResidue) = disordered_res[getdefaultresname(disordered_res)]
getresnames(disordered_res::DisorderedResidue) = sort(collect(keys(disordered_res.res_names)))

getresname(disordered_res::DisorderedResidue) = getdefaultresname(disordered_res)
getchainid(disordered_res::DisorderedResidue) = getchainid(getdefaultresidue(disordered_res))
getresnumber(disordered_res::DisorderedResidue) = getresnumber(getdefaultresidue(disordered_res))
getinscode(disordered_res::DisorderedResidue) = getinscode(getdefaultresidue(disordered_res))
ishetres(disordered_res::DisorderedResidue) = ishetres(getdefaultresidue(disordered_res))
geatomnames(disordered_res::DisorderedResidue) = getatomnames(getdefaultresidue(disordered_res))

getresid(disordered_res::DisorderedResidue) = getresid(getdefaultresidue(disordered_res))


getchainid(chain::Chain) = chain.id

function getresids(chain::Chain)

end


getmodelnumber(model::Model) = model.number
getchainids(model::Model) = sort(collect(keys(model.chains)))


getstrucname(struc::Structure) = struc.name
getmodelnumbers(struc::Structure) = sort(collect(keys(struc.models)))



start(::Structure) = 1

function next(struc::Structure, state)
    return (struc[getmodelnumbers(struc)[state]], state + 1)
end

done(struc::Structure, state) = state > length(getmodelnumbers(struc))

start(::Model) = 1

function next(model::Model, state)
    return (model[getchainids(model)[state]], state + 1)
end

done(model::Model, state) = state > length(getchainids(model))

start(::Chain) = 1

function next(chain::Chain, state)
    residues = sort(collect(values(chain.residues)), by= res -> res.ins_code)
    sort!(residues, by= res -> res.res_no)
    sort!(residues, by= res -> res.het_res)
    return (residues[state], state + 1)
end

done(chain::Chain, state) = state > length(getresids(chain))

start(::Residue) = 1

function next(res::Residue, state)
    return (res[getatomnames(res)[state]], state + 1)
end

done(res::Residue, state) = state > length(getatomnames(res))

start(::DisorderedResidue) = 1

function next(disordered_res::DisorderedResidue, state)
    reuturn (getdefaultresidue(disordered_res)[getatomnames(disordered_res)[state]], state + 1)
end

done(disordered_res::DisorderedResidue, state) = state > length(getatomnames(disordered_res))

start(::Atom) = 1

function next(atom::Atom, state)
    return (atom, state + 1)
end

done(atom::Atom, state) = state > 1

start(::DisorderedAtom) = 1

function next(disordered_atom::DisorderedAtom, state)
    return (disordered_atom[getaltlocs(disordered_atom)[state]], state + 1)
end

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

unfoldatoms(struc::Structure, args...) = [unfoldatoms(model, args...) for model in struc]
unfoldatoms(model::Model, args...) = [unfoldatoms(chain, args...) for chain in model]
unfoldatoms(chain::Chain, args...) = [unfoldatoms(res, args...) for res in chain]
unfoldatoms(res::AbstractResidue, args...) = unfold(res, args...)
unfoldatoms(atom::AbstractAtom, args...) = unfold(atom, args...)


# Do we need to check we are not overwriting? Yes
organise(models::ModelList; struc_name::AbstractString="") = Struc(struc_name, [getmodelnumber(model) => model for model in models])
organise(chains::ChainList; model_number::Int=1) = Model(model_number, [getchainid(chain) => chain for chain in chains])
organise(residues::ResidueList) =
organise(atoms::AtomList) =

organise(model::Model) = organise([model])
organise(chain::Chain) = organise([chain])
organise(res::AbstractResidue) = organise([res])
organise(atom::AbstractAtom, args...) = organise([atom], args...)

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
                    "CZ3", "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1",
                    "NH2", "NZ", "O", "OD1", "OD2", "OE1", "OE2", "OG",
                    "OG1", "OH", "SG"]

heavyatomselector(atom::AbstractAtom) = atomnameselector(atom, heavy_atom_names)

resnameselector(atom::AbstractAtom, res_names::Array{AbstractString,1}) = getresname(atom) in res_names

const water_res_names = []

waterselector(atom::AbstractAtom) = resnameselector(atom, water_res_names)
