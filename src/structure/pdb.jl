export PDB,
    PDBParseError,
    downloadpdb,
    parseatomrecord,
    spaceatomname,
    pdbline,
    writepdb,
    writepdblines


import Bio.FileFormat


"""Protein Data Bank (PDB) file format."""
immutable PDB <: FileFormat end


"""Error arising from parsing a PDB file."""
type PDBParseError <: Exception
    message::ASCIIString
    line_number::Int
    line::ASCIIString
end


Base.showerror(io::IO, e::PDBParseError) = println(io, e.message, " at line ", e.line_number, " of file:\n", e.line)


"""Download a PDB file or biological assembly from the RCSB PDB."""
function downloadpdb(pdbid::ASCIIString, out_filepath::ASCIIString="$pdbid.pdb"; ba_number::Int=0)
    # Check PDB ID is 4 characters long and only consits of alphanumeric characters
    @assert length(pdbid) == 4 && !ismatch(r"[^a-zA-Z0-9]", pdbid) "Not a valid PDB ID: \"$pdbid\""
    if ba_number == 0
        download("http://www.rcsb.org/pdb/files/$pdbid.pdb", out_filepath)
    else
        download("http://www.rcsb.org/pdb/files/$pdbid.pdb$ba_number", out_filepath)
    end
end


# Add selectors back in
function Base.read(input::IO,
            ::Type{PDB},
            selector_functions::Function...;
            structure_name::ASCIIString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true)
    # Dictionary of model numbers and raw atom lists
    atom_lists = Dict(1 => Atom[])
    # Entries outside of a MODEL/ENDMDL block are added to model 1
    curr_model = 1
    line_number = 0
    for line in eachline(input)
        line_number += 1
        # Read ATOM and HETATM records as required
        if (read_std_atoms && startswith(line, "ATOM  ")) || (read_het_atoms && startswith(line, "HETATM"))
            push!(atom_lists[curr_model], parseatomrecord(rstrip(line, '\n'), line_number))
        # Read MODEL record
        elseif startswith(line, "MODEL ")
            try
                curr_model = parse(Int, line[11:14])
            catch ArgumentError
                throw(PDBParseError("Could not read model serial number", line_number, line))
            end
            # Create model if required
            !haskey(atom_lists, curr_model) ? atom_lists[curr_model] = AbstractAtom[] : nothing
        # Read ENDMDL record
        elseif startswith(line, "ENDMDL")
            curr_model = 1
        end
    end
    # Remove model 1 atom list if it was not added to
    length(atom_lists[1]) == 0 ? delete!(atom_lists, 1) : nothing
    # Form disordered atom containers or remove atoms depending on remove_disorder
    # Apply selectors at this point (so e.g. disorderselector work correctly)
    # Form structure by organising each atom list into a model
    return organisestructure(
        [organisemodel(applyselectors(formatomlist(atom_list; remove_disorder=remove_disorder), selector_functions...); model_number=model_number) for (model_number, atom_list) in atom_lists]; structure_name=structure_name
    )
end


# Add selectors back in
function Base.read(filepath::ASCIIString,
            ::Type{PDB},
            selector_functions::Function...;
            structure_name::ASCIIString=splitdir(filepath)[2],
            kwargs...)
    open(filepath, "r") do input
        read(input, PDB, selector_functions...; structure_name=structure_name, kwargs...)
    end
end


"""Parse a PDB ATOM or HETATM record and return an `Atom`."""
function parseatomrecord(line::ASCIIString, line_number::Int=1)
    @assert startswith(line, "ATOM  ") || startswith(line, "HETATM") "Line does not appear to be an ATOM/HETATM record: \"$line\""
    return Atom(
        line[1:6] == "HETATM",
        parsestrict(line, (7,11), Int, "Could not read atom serial number", line_number),
        strip(parsestrict(line, (13,16), ASCIIString, "Could not read atom name", line_number)),
        parsestrict(line, (17,17), Char, "Could not read alt loc identifier", line_number),
        strip(parsestrict(line, (18,20), ASCIIString, "Could not read residue name", line_number)),
        parsestrict(line, (22,22), Char, "Could not read chain ID", line_number),
        parsestrict(line, (23,26), Int, "Could not read residue number", line_number),
        parsestrict(line, (27,27), Char, "Could not read insertion code", line_number),
        [
            parsestrict(line, (31,38), Float64, "Could not read x coordinate", line_number),
            parsestrict(line, (39,46), Float64, "Could not read y coordinate", line_number),
            parsestrict(line, (47,54), Float64, "Could not read z coordinate", line_number)
        ],
        parselenient(line, (55,60), Float64, 1.0),
        parselenient(line, (61,66), Float64, 0.0),
        strip(parselenient(line, (77,78), ASCIIString, "")),
        strip(parselenient(line, (79,80), ASCIIString, ""))
    )
end


"""Parse columns from a line and return the value or throw an error."""
function parsestrict(line::ASCIIString,
                    cols::Tuple{Int, Int},
                    out_type::Type,
                    error_message::ASCIIString,
                    line_number::Int)
    try
        return parsevalue(line, cols, out_type)
    catch
        throw(PDBParseError(error_message, line_number, line))
    end
end


"""Parse columns from a line and return the value or a default value."""
function parselenient(line::ASCIIString,
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
function parsevalue(line::ASCIIString, cols::Tuple{Int, Int}, out_type::Type)
    try
        if out_type == Int
            return parse(Int, line[cols[1]:cols[2]])
        elseif out_type == Float64
            return parse(Float64, line[cols[1]:cols[2]])
        elseif out_type == ASCIIString
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
        return string_out[1:new_length]
    else
        return "$(repeat(" ", new_length - length(string_out)))$string_out"
    end
end


"""Space an `Atom` name such that the last element letter (generally) appears in the second
column. If the `element` property of the `Atom` is set it is used to get the element,
otherwise the name starts from the second column where possible."""
function spaceatomname(atom::Atom)
    atom_name = atomname(atom)
    chars = length(atom_name)
    @assert chars <= 4 "Atom name is greater than four characters: \"$atom_name\""
    # In the absence of the element, the first index goes in column two
    if element(atom) == "" || findfirst(atom_name, element(atom)[1]) == 0
        cent_ind = 1
    # The last letter of the element goes in column two where possible
    else
        cent_ind = findfirst(atom_name, element(atom)[1]) + length(element(atom)) - 1
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


"""Form a Protein Data Bank (PDB) format ATOM/HETATM record from an `Atom`."""
pdbline(atom::Atom) = ASCIIString[
        ishetatom(atom) ? "HETATM" : "ATOM  ",
        spacestring(serial(atom), 5),
        " ",
        spaceatomname(atom),
        string(altlocid(atom)),
        spacestring(resname(atom), 3),
        " ",
        string(chainid(atom)),
        spacestring(resnumber(atom), 4),
        string(inscode(atom)),
        "   ",
        spacestring(round(x(atom), 3), 8),
        spacestring(round(y(atom), 3), 8),
        spacestring(round(z(atom), 3), 8),
        "  ",
        spacestring(round(occupancy(atom), 2), 4),
        " ",
        spacestring(round(tempfac(atom), 2), 5),
        "          ",
        spacestring(element(atom), 2),
        spacestring(charge(atom), 2),
    ]


"""Write a `StructuralElementOrList` to a PDB format file."""
function writepdb(output::IO, element::Union{ProteinStructure, Vector{Model}}, selector_functions::Function...)
    # If there are multiple models, write out MODEL/ENDMDL lines
    if length(element) > 1
        for model in element
            println(output, "MODEL     ", spacestring(modelnumber(model), 4), repeat(" ", 66))
            writepdblines(output, model, selector_functions...)
            println(output, "ENDMDL$(repeat(" ", 74))")
        end
    # If there is only one model, do not write out MODEL/ENDMDL lines
    else
        writepdblines(output, element, selector_functions...)
    end
end

writepdb(output::IO, element::StructuralElementOrList, selector_functions::Function...) = writepdblines(output, element, selector_functions...)


"""Write a `StructuralElementOrList` as lines in PDB format."""
function writepdblines(output::IO, element::StructuralElementOrList, selector_functions::Function...)
    for atom in collectatoms(element, selector_functions...), atom_record in atom
        println(output, pdbline(atom_record)...)
    end
end


function writepdb(filepath::ASCIIString, element::StructuralElementOrList, selector_functions::Function...)
    open(filepath, "w") do output
        writepdb(output, element, selector_functions...)
    end
end
