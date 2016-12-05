export
    PDB,
    PDBParseError,
    downloadpdb,
    spaceatomname,
    pdbline,
    writepdb


"Protein Data Bank (PDB) file format."
immutable PDB <: FileFormat end


"Error arising from parsing a Protein Data Bank (PDB) file."
type PDBParseError <: Exception
    message::String
    line_number::Int
    line::String
end


function Base.showerror(io::IO, e::PDBParseError)
    return print(io,
            e.message,
            " at line ",
            e.line_number,
            " of file:\n",
            e.line)
end


"""
Download a Protein Data Bank (PDB) file or biological assembly from the RCSB
PDB. By default downloads the PDB file; if `ba_number` is set the biological
assembly with that number will be downloaded.
"""
function downloadpdb(pdbid::AbstractString,
                    out_filepath::AbstractString="$pdbid.pdb";
                    ba_number::Integer=0)
    # Check PDB ID is 4 characters long and only consits of alphanumeric characters
    if length(pdbid) != 4 || ismatch(r"[^a-zA-Z0-9]", pdbid)
        throw(ArgumentError("Not a valid PDB ID: \"$pdbid\""))
    end
    if ba_number == 0
        download("http://www.rcsb.org/pdb/files/$pdbid.pdb", out_filepath)
    else
        # Will download error page if ba_number is too high
        download("http://www.rcsb.org/pdb/files/$pdbid.pdb$ba_number", out_filepath)
    end
end


function Base.read(input::IO,
            ::Type{PDB};
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true)
    # Define ProteinStructure and add to it incrementally
    struc = ProteinStructure(structure_name)
    struc[1] = Model(1, struc)
    # Entries outside of a MODEL/ENDMDL block are added to model 1
    curr_model = 1
    line_n = 0
    for line in eachline(input)
        line_n += 1
        # Read ATOM and HETATM records as required
        if (read_std_atoms && startswith(line, "ATOM  ")) ||
                (read_het_atoms && startswith(line, "HETATM"))
            unsafe_addatomtomodel!(
                    struc[curr_model],
                    AtomRecord(line, line_n), remove_disorder=remove_disorder)
        # Read MODEL record
        elseif startswith(line, "MODEL ")
            try
                curr_model = parse(Int, line[11:14])
            catch
                throw(PDBParseError(
                    "Could not read model serial number", line_n, line))
            end
            # Create model if required
            if !haskey(models(struc), curr_model)
                struc[curr_model] = Model(curr_model, struc)
            end
        # Read ENDMDL record
        elseif startswith(line, "ENDMDL")
            curr_model = 1
        end
    end
    # Remove any models that were not added to
    for model in struc
        if countchains(model) == 0
            delete!(models(struc), modelnumber(model))
        end
    end
    # Generate lists for iteration
    fixlists!(struc)
    return struc
end

function Base.read(filepath::AbstractString,
            ::Type{PDB};
            structure_name::AbstractString=splitdir(filepath)[2],
            kwargs...)
    open(filepath, "r") do input
        read(input, PDB; structure_name=structure_name, kwargs...)
    end
end


# Constructor from PDB ATOM/HETATM line
AtomRecord(pdb_line::String, line_n::Integer=1) = AtomRecord(
    pdb_line[1] == 'H', # This assumes the line has already been checked as an ATOM/HETATM record
    parsestrict(pdb_line, 7, 11, Int, "Could not read atom serial number", line_n),
    parsestrict(pdb_line, 13, 16, String, "Could not read atom name", line_n),
    parsestrict(pdb_line, 17, 17, Char, "Could not read alt loc identifier", line_n),
    parsestrict(pdb_line, 18, 20, String, "Could not read residue name", line_n),
    parsestrict(pdb_line, 22, 22, Char, "Could not read chain ID", line_n),
    parsestrict(pdb_line, 23, 26, Int, "Could not read residue number", line_n),
    parsestrict(pdb_line, 27, 27, Char, "Could not read insertion code", line_n),
    [
        parsestrict(pdb_line, 31, 38, Float64, "Could not read x coordinate", line_n),
        parsestrict(pdb_line, 39, 46, Float64, "Could not read y coordinate", line_n),
        parsestrict(pdb_line, 47, 54, Float64, "Could not read z coordinate", line_n)
    ],
    parselenient(pdb_line, 55, 60, Float64, 1.0),
    parselenient(pdb_line, 61, 66, Float64, 0.0),
    parselenient(pdb_line, 77, 78, String, "  "),
    parselenient(pdb_line, 79, 80, String, "  ")
)


"Parse columns from a line and return the value or throw a `PDBParseError`."
function parsestrict(line::String,
                    col_1::Integer,
                    col_2::Integer,
                    out_type::Type,
                    error_message::AbstractString,
                    line_n::Integer)
    try
        return parsevalue(line, col_1, col_2, out_type)
    catch
        throw(PDBParseError(error_message, line_n, line))
    end
end


"Parse columns from a line and return the value or a default value."
function parselenient(line::String,
                    col_1::Integer,
                    col_2::Integer,
                    out_type::Type,
                    default)
    try
        return parsevalue(line, col_1, col_2, out_type)
    catch
        return default
    end
end


"Parse columns from a line into a given type."
function parsevalue(line::String,
                    col_1::Integer,
                    col_2::Integer,
                    out_type::Type)
    try
        if out_type == Int
            return parse(Int, line[col_1:col_2])
        elseif out_type == Float64
            return parse(Float64, line[col_1:col_2])
        elseif out_type == String
            return line[col_1:col_2]
        elseif out_type == Char
            return line[col_1]
        else
            error()
        end
    catch
        error("Could not parse to desired type")
    end
end


"""
Form a string of a certain length from a value by adding spaces to the left.
Throws an error if the value is too long.
"""
function spacestring(val_in, new_length::Integer)
    string_out = string(val_in)
    if length(string_out) > new_length
        throw(ArgumentError("Cannot fit value \"$string_out\" into $new_length space(s)"))
    end
    return lpad(string_out, new_length)
end


"""
Space an `Atom` name such that the last element letter (generally) appears in
the second column. If the `element` property of the `Atom` is set it is used to
get the element, otherwise the name starts from the second column where
possible. This function is generally not required as spacing is recorded when
atom names are read in from a Protein Data Bank (PDB) file.
"""
function spaceatomname(at::Atom)
    at_name = atomname(at, strip=false)
    chars = length(at_name)
    if chars == 4
        return at_name
    end
    strip_el = element(at, strip=true)
    if chars > 4
        throw(ArgumentError("Atom name is greater than four characters: \"$at_name\""))
    end
    # In the absence of the element, the first index goes in column two
    if strip_el == "" || findfirst(at_name, strip_el[1]) == 0
        cent_ind = 1
    # The last letter of the element goes in column two where possible
    else
        cent_ind = findfirst(at_name, strip_el[1]) + length(strip_el) - 1
    end
    if cent_ind > 2
        throw(ArgumentError("Atom name is too long to space correctly: \"$at_name\""))
    end
    if cent_ind == 1 && chars < 4
        out_string = " $at_name"
    else
        out_string = "$at_name"
    end
    return rpad(out_string, 4)
end


"""
Form a Protein Data Bank (PDB) format ATOM or HETATM record from an `Atom` or
`AtomRecord`.
"""
function pdbline(at::Atom)
    return (ishetero(at) ? "HETATM" : "ATOM  ") *
            spacestring(serial(at), 5) *
            " " *
            spaceatomname(at) *
            string(altlocid(at)) *
            spacestring(resname(at), 3) *
            " " *
            string(chainid(at)) *
            spacestring(resnumber(at), 4) *
            string(inscode(at)) *
            "   " *
            # This will throw an error for large coordinate values, e.g. -1000.123
            spacestring(round(x(at), 3), 8) *
            spacestring(round(y(at), 3), 8) *
            spacestring(round(z(at), 3), 8) *
            spacestring(round(occupancy(at), 2), 6) *
            # This will throw an error for large temp facs, e.g. 1000.12
            spacestring(round(tempfactor(at), 2), 6) *
            "          " *
            spacestring(element(at), 2) *
            spacestring(charge(at), 2)
end

function pdbline(at_rec::AtomRecord)
    return (at_rec.het_atom ? "HETATM" : "ATOM  ") *
            spacestring(at_rec.serial, 5) *
            " " *
            spacestring(at_rec.atom_name, 4) *
            string(at_rec.alt_loc_id) *
            spacestring(at_rec.res_name, 3) *
            " " *
            string(at_rec.chain_id) *
            spacestring(at_rec.res_number, 4) *
            string(at_rec.ins_code) *
            "   " *
            # This will throw an error for large coordinate values, e.g. -1000.123
            spacestring(round(at_rec.coords[1], 3), 8) *
            spacestring(round(at_rec.coords[2], 3), 8) *
            spacestring(round(at_rec.coords[3], 3), 8) *
            spacestring(round(at_rec.occupancy, 2), 6) *
            # This will throw an error for large temp facs, e.g. 1000.12
            spacestring(round(at_rec.temp_factor, 2), 6) *
            "          " *
            spacestring(at_rec.element, 2) *
            spacestring(at_rec.charge, 2)
end


"""
Write a `StructuralElementOrList` to a Protein Data Bank (PDB) format file. Only
ATOM, HETATM, MODEL and ENDMDL records are written - there is no header and
there are no TER records.
Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function writepdb(output::IO,
                el::Union{ProteinStructure, Vector{Model}},
                atom_selectors::Function...)
    # If there are multiple models, write out MODEL/ENDMDL lines
    if length(el) > 1
        for mod in sort(collect(el))
            println(output, "MODEL     ", spacestring(modelnumber(mod), 4), repeat(" ", 66))
            writepdb(output, mod, atom_selectors...)
            println(output, "ENDMDL$(repeat(" ", 74))")
        end
    # If there is only one model, do not write out MODEL/ENDMDL lines
    elseif isa(el, ProteinStructure)
        writepdb(output, defaultmodel(el), atom_selectors...)
    else
        writepdb(output, el[1], atom_selectors...)
    end
end


function writepdb(output::IO,
                el::Union{Model, Chain, AbstractResidue, Vector{Chain},
                    Vector{AbstractResidue}, Vector{Residue}, Vector{DisorderedResidue}},
                atom_selectors::Function...)
    # Collect residues then expand out disordered residues
    for res in collectresidues(el)
        if isa(res, Residue)
            for at in collectatoms(res, atom_selectors...)
                writepdb(output, at)
            end
        else
            for res_name in resnames(res)
                for at in collectatoms(disorderedres(res, res_name), atom_selectors...)
                    writepdb(output, at)
                end
            end
        end
    end
end

function writepdb(output::IO, at::AbstractAtom, atom_selectors::Function...)
    for atom_record in at
        println(output, pdbline(atom_record))
    end
end

function writepdb{T <: AbstractAtom}(output::IO,
                    ats::Vector{T},
                    atom_selectors::Function...)
    for at in collectatoms(ats, atom_selectors...)
        writepdb(output, at)
    end
end

function writepdb(filepath::AbstractString,
                el::StructuralElementOrList,
                atom_selectors::Function...)
    open(filepath, "w") do output
        writepdb(output, el, atom_selectors...)
    end
end
