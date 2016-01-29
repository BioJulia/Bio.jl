import Base: showerror,
    read

import Bio.FileFormat


immutable PDB <: FileFormat end


type PDBParseException <: Exception
    message::AbstractString
    line_no::Int
    line::AbstractString
end


function showerror(io::IO, e::PDBParseException)
    println(io, e.message, " at line ", e.line_no, " of file:\n", e.line)
end


function getpdb(pdbid::AbstractString, out_filepath::AbstractString="$pdbid.pdb"; ba_number::Int=0)
    # Check PDB ID is 4 characters and only consits of alphanumeric characters
    @assert length(pdbid) == 4 && !ismatch(r"[^a-zA-Z0-9]", pdbid) "Not a valid PDB ID: \"$pdbid\""
    if ba_number == 0
        download("http://www.rcsb.org/pdb/files/$pdbid.pdb", out_filepath)
    else
        download("http://www.rcsb.org/pdb/files/$pdbid.pdb$ba_number", out_filepath)
    end
end


function read(input::IO,
            ::Type{PDB},
            args...;
            struc_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true)
    # Dictionary of model numbers and atom lists
    atom_lists = Dict(1 => AbstractAtom[])
    # Entries outside of a MODEL/ENDMDL block are added to model 1
    curr_model = 1
    line_no = 0
    for line in eachline(input)
        line_no += 1
        # Read ATOM and HETATM records as required
        if (read_std_atoms && startswith(line, "ATOM  ")) || (read_het_atoms && startswith(line, "HETATM"))
            atom = parseatomrecord(rstrip(line, '\n'), line_no)
            # Check selector functions are satisfied and if not skip this atom
            selected = true
            for selector_function in args
                !selector_function(atom) ? (selected = false; break)
            end
            !selected ? continue
            try
                addtoatomlist!(atom_lists[curr_model], atom, remove_disorder=remove_disorder)
            catch ErrorException
                throw(PDBParseException("Two copies of the same atom have the same alternative location ID (note names are stripped of whitespace so identical names with different spacing, e.g. ' CA ' and 'CA  ', are not currently supported)", line_no, line))
            end
        # Read MODEL record
        elseif startswith(line, "MODEL ")
            try
                curr_model = parse(Int, line[11:14])
            catch ArgumentError
                throw(PDBParseException("Could not read model serial number", line_no, line))
            end
            # Create model if required
            !haskey(atom_lists, curr_model) ? atom_lists[curr_model] = AbstractAtom[]
        # Read ENDMDL record
        elseif startswith(line, "ENDMDL")
            curr_model = 1
        end
    end
    # Remove model 1 atom list if it was not added to
    length(atom_lists[1]) == 0 ? delete!(atom_lists[1], 1)
    # Form structure by organising each atom list into a model
    return organisestruc([organisemodel(atom_list) for atom_list in values(atom_lists)], struc_name)
end


function read(filepath::AbstractString,
            ::Type{PDB},
            args...;
            struc_name::AbstractString=splitdir(filepath)[2],
            kwargs...)
    input = open(filepath)
    finalizer(input, input -> close(input))
    return read(input, PDB, args...; struc_name=struc_name, kwargs...)
end


function parseatomrecord(line::AbstractString, line_no::Int=1)
    @assert startswith(line, "ATOM  ") || startswith(line, "HETATM") "Line does not appear to be an ATOM/HETATM record: \"$line\""
    het_atom = line[1:6] == "HETATM"
    serial = parsestrict(line, (7,11), Int, "Could not read atom serial number", line_no)
    name = strip(parsestrict(line, (13,16), AbstractString, "Could not read atom name", line_no))
    alt_loc = parsestrict(line, (17,17), Char, "Could not read alt loc identifier", line_no)
    res_name = strip(parsestrict(line, (18,20), AbstractString, "Could not read residue name", line_no))
    chain_id = parsestrict(line, (22,22), Char, "Could not read chain ID", line_no)
    res_no = parsestrict(line, (23,26), Int, "Could not read residue number", line_no)
    ins_code = parsestrict(line, (27,27), Char, "Could not read insertion code", line_no)
    x = parsestrict(line, (31,38), Float64, "Could not read x coordinate", line_no)
    y = parsestrict(line, (39,46), Float64, "Could not read y coordinate", line_no)
    z = parsestrict(line, (47,54), Float64, "Could not read z coordinate", line_no)
    occupancy = parselenient(line, (55,60), Float64, 1.0)
    temp_fac = parselenient(line, (61,66), Float64, 0.0)
    element = strip(parselenient(line, (77,78), AbstractString, ""))
    charge = strip(parselenient(line, (79,80), AbstractString, ""))
    return Atom(
        het_atom, serial, name, alt_loc, res_name, chain_id, res_no, ins_code,
        [x, y, z], occupancy, temp_fac, element, charge
    )
end


"""Parse columns from a line and return the value or throw an error."""
function parsestrict(line::AbstractString,
                    cols::Tuple{Int, Int},
                    out_type::Type,
                    error_message::AbstractString,
                    line_no::Int)
    try
        return parsevalue(line, cols, out_type)
    catch
        throw(PDBParseException(error_message, line_no, line))
    end
end


"""Parse columns from a line and return the value or a default value."""
function parselenient(line::AbstractString,
                    cols::Tuple{Int, Int},
                    out_type::Type,
                    default)
    try
        return parsevalue(line, cols, out_type)
    catch
        return default
    end
end


"""Parse columns from a line."""
function parsevalue(line::AbstractString, cols::Tuple{Int, Int}, out_type::Type)
    try
        if out_type == Int
            return parse(Int, line[cols[1]:cols[2]])
        elseif out_type == Float64
            return parse(Float64, line[cols[1]:cols[2]])
        elseif out_type == AbstractString
            return line[cols[1]:cols[2]]
        elseif out_type == Char
            return line[cols[1]]
        else
            error()
        end
    catch
        error("Could not parse to desired type")
    end
end


"""Form a string of a certain length from a value."""
function spacestring(val_in, new_length::Int)
    string_out = string(val_in)
    if length(string_out) > new_length
        string_out = string_out[1:new_length]
    else
        string_out = "$(repeat(" ", new_length - length(string_out)))$string_out"
    end
    return string_out
end


"""
Space an `Atom` name such that the second element letter (generally) appears in the second
column. Having the `element` property of the `Atom` set improves the result.
"""
function spaceatomname(atom::Atom)
    atom_name = atom.name
    chars = length(atom_name)
    @assert chars <= 4 "Atom name is greater than four characters: \"$atom_name\""
    if atom.element == "" || findfirst(atom_name, atom.element[1]) == 0
        # In the absence of the element, the first index goes in column two
        cent_ind = 1
    else
        # The last letter of the element goes in column two where possible
        cent_ind = findfirst(atom_name, atom.element[1]) + length(atom.element) - 1
    end
    @assert cent_ind <= 2 "Atom name is too long to space correctly: \"$atom_name\""
    if cent_ind == 1 && chars < 4
        out_string = " $atom_name"
    else
        out_string = "$atom_name"
    end
    while length(out_string) < 4
        out_string = "$out_string "
    end
    return out_string
end


"""
Form a Protein Data Bank (PDB) format ATOM/HETATM record from an `Atom`.
"""
getpdbline(atom::Atom) = ASCIIString[
        ishetatom(atom) ? "HETATM" : "ATOM  ",
        spacestring(getserial(atom), 5),
        " ",
        spaceatomname(atom),
        string(getaltlocid(atom)),
        spacestring(getresname(atom), 3),
        " ",
        string(getchainid(atom)),
        spacestring(getresnumber(atom), 4),
        string(getinscode(atom)),
        "   ",
        spacestring(round(getx(atom), 3), 8),
        spacestring(round(gety(atom), 3), 8),
        spacestring(round(getz(atom), 3), 8),
        "  ",
        spacestring(round(getoccupancy(atom), 2), 4),
        " ",
        spacestring(round(gettempfac(atom), 2), 5),
        "          ",
        spacestring(getelement(atom), 2),
        spacestring(getcharge(atom), 2),
    ]


function writepdb(output::IO, element::Union{Structure, ModelList}, args...)
    # If there are multiple models, write out MODEL/ENDMDL lines
    if length(element) > 1
        for model in element
            println(output, "MODEL     ", spacestring(getmodelnumber(model), 4), repeat(" ", 66))
            writepdblines(output, model)
            println(output, "ENDMDL$(repeat(" ", 74))")
        end
    # If there is only one model, do not write out MODEL/ENDMDL lines
    else
        writepdblines(output, element, args...)
    end
end

writepdb(output::IO, element::StrucElementOrList, args...) = writepdblines(output, element, args...)


function writepdblines(output::IO, element::StrucElementOrList, args...)
    for atom in collectatoms(element, args...), atom_record in atom
        println(output, getpdbline(atom_record)...)
    end
end


function writepdb(filepath::AbstractString, element::StrucElementOrList, args...)
    output = open(filepath, "w")
    finalizer(output, output -> close(output))
    return writepdb(output, element, args...)
end
