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

typealias AtomList Array{AbstractAtom,1}

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

typealias ResidueList Array{AbstractResidue,1}

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


typealias StrucElementOrList Union{StrucElement, AtomList, ResidueList}


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



start(::Structure) = 1

function next(struc::Structure, state)
    return (struc.models[sort(collect(keys(struc.models)))[state]], state + 1)
end

done(struc::Structure, state) = state > length(struc.models)

start(::Model) = 1

function next(model::Model, state)
    return (model.chains[sort(collect(keys(model.chains)))[state]], state + 1)
end

done(model::Model, state) = state > length(model.chains)

start(::Chain) = 1

function next(chain::Chain, state)
    residues = sort(collect(values(chain.residues)), by= res -> res.ins_code)
    sort!(residues, by= res -> res.res_no)
    sort!(residues, by= res -> res.het_res)
    return (residues[state], state + 1)
end

done(chain::Chain, state) = state > length(chain.residues)

start(::Residue) = 1

function next(res::Residue, state)
    return (sort(collect(values(res.atoms)), by= atom -> atom.serial)[state], state + 1)
end

done(res::Residue, state) = state > length(res.atoms)

start(::DisorderedResidue) = 1

function next(disordered_res::DisorderedResidue, state)
    return (sort(collect(values(disordered_res[disordered_res.default].atoms)), by= atom -> atom.serial)[state], state + 1)
end

done(disordered_res::DisorderedResidue, state) = state > length(disordered_res[disordered_res.default].atoms)

start(::Atom) = 1

function next(atom::Atom, state)
    return (atom, state + 1)
end

done(atom::Atom, state) = state > 1

start(::DisorderedAtom) = 1

function next(disordered_atom::DisorderedAtom, state)
    return (disordered_atom.alt_locs[sort(collect(keys(disordered_atom.alt_locs)))[state]], state + 1)
end

done(disordered_atom::DisorderedAtom, state) = state > length(disordered_atom.alt_locs)


function filteratoms!(atoms::AtomList, args...)
    for selection_function in args
        filter!(atoms, atom -> selection_function(atom))
    end
end


function unfold(element::StrucElement)
    return [sub_element for sub_element in element]
end

function unfold(atom_container::Union{AbstractResidue, AbstractAtom}, args...)
    atoms = unfold(atom_container)
    filteratoms!(atoms, args...)
    return atoms
end

function unfoldresidues()

end

function unfoldatoms(, args...)

end


function organise()

end

function organisestruc()

end

function organisemodel()

end


stdatomselector(atom::AbstractAtom) = !ishetatom(atom)

hetatomselector(atom::AbstractAtom) = ishetatom(atom)

function atomnameselector(atom::AbstractAtom, atom_name_list::Array{AbstractString,1})
    if getatomname(atom) in atom_name_list
        return true
    else
        return false
    end
end

const calpha_list = ["CA"]

calphaselector(atom::AbstractAtom) = atomnameselector(atom, calpha_list)

const backbone_list = ["CA", "N", "C"]

backboneselector(atom::AbstractAtom) = atomnameselector(atom, backbone_list)

const heavy_atom_list = ["C", "CA", "CB", "CD", "CD1", "CD2", "CE", "CE1",
                    "CE2", "CE3", "CG", "CG1", "CG2", "CH2", "CZ", "CZ2",
                    "CZ3", "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1",
                    "NH2", "NZ", "O", "OD1", "OD2", "OE1", "OE2", "OG",
                    "OG1", "OH", "SG"]

heavyatomselector(atom::AbstractAtom) = atomnameselector(atom, heavy_atom_list)
