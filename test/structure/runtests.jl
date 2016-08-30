module TestStructure

using Base.Test

using Bio.Structure
using TestFunctions.get_bio_fmt_specimens
using Bio.Structure:
    atomid,
    fixlists!,
    parsestrict,
    parselenient,
    parsevalue,
    spacestring,
    writepdblines


get_bio_fmt_specimens()

# Access a PDB file in BioFmtSpecimens
pdbfilepath(filename::AbstractString) = Pkg.dir("Bio", "test", "BioFmtSpecimens", "PDB", filename)


@testset "Model" begin
    # Test constructors and indexing
    struc = ProteinStructure("Test structure")
    struc[1] = Model(1, struc)
    mod = struc[1]
    @test isa(mod, Model)
    struc[3] = Model(3, struc)
    struc[1]['A'] = Chain('A', mod)
    ch = struc[1]['A']
    @test isa(ch, Chain)
    struc['B'] = Chain('B', mod)
    struc['A'][10] = Residue("ALA", 10, ' ', false, ch)
    res = struc['A'][10]
    @test isa(res, Residue)
    struc['A']["H_20A"] = DisorderedResidue(Dict(
        "VAL"=> Residue("VAL", 20, 'A', true, ch),
        "ILE"=> Residue("ILE", 20, 'A', true, ch)
    ), "VAL")
    dis_res = struc['A']["H_20A"]
    @test isa(dis_res, DisorderedResidue)
    struc['A'][10][" CA "] = Atom(100, " CA ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res)
    at = struc['A'][10]["CA"]
    @test isa(at, Atom)
    struc['A'][10][" CB "] = DisorderedAtom(Dict(
        'A'=> Atom(200, " CB ", 'A', [10.0, 20.0, 30.0], 0.6, 20.0, " C", "  ", res),
        'B'=> Atom(201, " CB ", 'B', [11.0, 21.0, 31.0], 0.4, 30.0, " C", "  ", res)
    ), 'A')
    dis_at = struc['A'][10]["CB"]
    @test isa(dis_at, DisorderedAtom)
    struc['A']["H_20A"][" CG "] = Atom(300, " CG ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", defaultresidue(dis_res))
    disorderedres(dis_res, "ILE")[" O  "] = Atom(400, " O  ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " O", "  ", disorderedres(dis_res, "ILE"))
    fixlists!(struc)

    # Test show
    show(DevNull, at)
    showcompact(DevNull, at)
    show(DevNull, dis_at)
    show(DevNull, res)
    show(DevNull, dis_res)
    show(DevNull, ch)
    show(DevNull, mod)
    show(DevNull, struc)

    # Test getters/setters
    @test !ishetatom(at)
    @test !ishetatom(dis_at)

    @test serial(at) == 100
    @test serial(dis_at) == 200

    @test atomname(at) == "CA"
    @test atomname(dis_at) == "CB"
    @test atomname(at, spaces=true) == " CA "
    @test atomname(dis_at, spaces=true) == " CB "

    @test altlocid(at) == ' '
    @test altlocid(dis_at) == 'A'

    @test resname(at) == "ALA"
    @test resname(dis_at) == "ALA"

    @test chainid(at) == 'A'
    @test chainid(dis_at) == 'A'

    @test resnumber(at) == 10
    @test resnumber(dis_at) == 10

    @test inscode(at) == ' '
    @test inscode(dis_at) == ' '

    @test x(at) == 1.0
    @test x(dis_at) == 10.0

    @test y(at) == 2.0
    @test y(dis_at) == 20.0

    @test z(at) == 3.0
    @test z(dis_at) == 30.0

    @test coords(at) == [1.0, 2.0, 3.0]
    @test coords(dis_at) == [10.0, 20.0, 30.0]

    @test occupancy(at) == 1.0
    @test occupancy(dis_at) == 0.6

    @test tempfac(at) == 10.0
    @test tempfac(dis_at) == 20.0

    @test element(at) == "C"
    @test element(dis_at) == "C"
    @test element(at, spaces=true) == " C"
    @test element(dis_at, spaces=true) == " C"

    @test charge(at) == "  "
    @test charge(dis_at) == "  "

    @test !ishetero(at)
    @test !ishetero(dis_at)

    @test !isdisorderedatom(at)
    @test isdisorderedatom(dis_at)

    @test atomid(at) == ("10:A", "ALA", "CA")
    @test atomid(dis_at) == ("10:A", "ALA", "CB")

    @test resid(at) == "10"
    @test resid(dis_at) == "10"
    @test resid(at, full=true) == "10:A"
    @test resid(dis_at, full=true) == "10:A"

    x!(at, 5.0)
    @test coords(at) == [5.0, 2.0, 3.0]
    y!(at, 6.0)
    @test coords(at) == [5.0, 6.0, 3.0]
    z!(at, 7.0)
    @test coords(at) == [5.0, 6.0, 7.0]
    coords!(at, [40.0, 50.0, 60.0])
    @test coords(at) == [40.0, 50.0, 60.0]

    # Only the coordinates on the default atom are changed
    x!(dis_at, 12.0)
    @test coords(dis_at) == [12.0, 20.0, 30.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]
    y!(dis_at, 22.0)
    @test coords(dis_at) == [12.0, 22.0, 30.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]
    z!(dis_at, 32.0)
    @test coords(dis_at) == [12.0, 22.0, 32.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]
    coords!(dis_at, [40.0, 50.0, 60.0])
    @test coords(dis_at) == [40.0, 50.0, 60.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]
    x!(dis_at['A'], 100.0)
    @test coords(dis_at) == [100.0, 50.0, 60.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]
    x!(dis_at['B'], 110.0)
    @test coords(dis_at) == [100.0, 50.0, 60.0]
    @test coords(dis_at['B']) == [110.0, 21.0, 31.0]

    # Test Atom iteration
    atom_list = collect(at)
    @test isa(atom_list, Vector{Atom})
    @test length(atom_list) == 1
    @test serial(atom_list[1]) == 100

    # Test DisorderedAtom getters/setters
    @test defaultaltlocid(dis_at) == 'A'
    @test isa(defaultatom(dis_at), Atom)
    @test serial(defaultatom(dis_at)) == 200
    @test altlocids(dis_at) == ['A', 'B']

    dis_at_mod = DisorderedAtom(dis_at, 'B')
    @test defaultaltlocid(dis_at_mod) == 'B'
    @test serial(dis_at_mod) == 201
    @test_throws AssertionError DisorderedAtom(dis_at, 'C')

    # Test DisorderedAtom indices
    @test isa(dis_at['A'], Atom)
    @test serial(dis_at['A']) == 200
    @test serial(dis_at['B']) == 201
    @test_throws KeyError dis_at['C']
    @test_throws MethodError dis_at["A"]

    # Test DisorderedAtom iteration
    dis_at_list = collect(dis_at)
    @test isa(dis_at_list, Vector{Atom})
    @test length(dis_at_list) == 2
    @test serial(dis_at_list[2]) == 201

    @test resname(res) == "ALA"
    @test chainid(res) == 'A'
    @test resnumber(res) == 10
    @test inscode(res) == ' '
    @test !ishetres(res)
    @test atomnames(res) == ["CA", "CB"]
    @test isa(atoms(res), Dict{String, AbstractAtom})
    @test length(atoms(res)) == 2
    @test serial(atoms(res)[" CA "]) == 100

    @test !ishetero(res)
    @test !isdisorderedres(res)
    @test resid(res) == "10"
    @test resid(res, full=true) == "10:A"

    # Test Residue indices
    @test isa(res["CA"], AbstractAtom)
    @test serial(res["CA"]) == 100
    @test serial(res["CB"]) == 200
    @test_throws KeyError res["N"]

    # Test Residue iteration
    res_list = collect(res)
    @test isa(res_list, Vector{AbstractAtom})
    @test length(res_list) == 2
    @test serial(res_list[2]) == 200

    # Test DisorderedResidue getters/setters
    @test isa(disorderedres(dis_res, "ILE"), AbstractResidue)
    @test resname(disorderedres(dis_res, "ILE")) == "ILE"
    @test defaultresname(dis_res) == "VAL"
    @test isa(defaultresidue(dis_res), Residue)
    @test resname(defaultresidue(dis_res)) == "VAL"
    @test resnames(dis_res) == ["VAL", "ILE"]

    @test resname(dis_res) == "VAL"
    @test chainid(dis_res) == 'A'
    @test resnumber(dis_res) == 20
    @test inscode(dis_res) == 'A'
    @test ishetres(dis_res)

    @test atomnames(dis_res) == ["CG"]
    @test isa(atoms(dis_res), Dict{String, AbstractAtom})
    @test length(atoms(dis_res)) == 1
    @test serial(atoms(dis_res)[" CG "]) == 300

    @test ishetero(dis_res)
    @test isdisorderedres(dis_res)
    @test resid(dis_res) == "H_20A"
    @test resid(dis_res, full=true) == "H_20A:A"

    dis_res_mod = DisorderedResidue(dis_res, "ILE")
    @test defaultresname(dis_res_mod) == "ILE"
    @test atomnames(dis_res_mod) == ["O"]
    @test_throws AssertionError DisorderedResidue(dis_res, "SER")

    # Test DisorderedResidue indices
    @test isa(dis_res["CG"], AbstractAtom)
    @test serial(dis_res["CG"]) == 300
    @test_throws KeyError dis_res["N"]

    # Test DisorderedResidue iteration
    dis_res_list = collect(dis_res)
    @test isa(dis_res_list, Vector{AbstractAtom})
    @test length(dis_res_list) == 1
    @test serial(dis_res_list[1]) == 300

    # Test Chain getters/setters
    @test chainid(ch) == 'A'
    @test resids(ch) == ["10", "H_20A"]
    @test isa(residues(ch), Dict{String, AbstractResidue})
    @test length(residues(ch)) == 2
    @test serial(residues(ch)["10"]["CA"]) == 100

    # Test Chain indices
    @test isa(ch["10"], AbstractResidue)
    @test isa(ch[10], AbstractResidue)
    @test_throws KeyError ch["H_10"]
    @test serial(ch["10"]["CA"]) == 100
    @test serial(ch[10]["CA"]) == 100
    @test resname(ch["H_20A"]) == "VAL"
    @test_throws KeyError ch["11"]
    @test_throws KeyError ch[11]

    # Test Chain iteration
    chain_list = collect(ch)
    @test isa(chain_list, Vector{AbstractResidue})
    @test length(chain_list) == 2
    @test resname(chain_list[2]) == "VAL"

    # Test Model getters/setters
    @test modelnumber(mod) == 1
    @test chainids(mod) == ['A', 'B']
    @test isa(chains(mod), Dict{Char, Chain})
    @test length(chains(mod)) == 2
    @test resname(chains(mod)['A']["H_20A"]) == "VAL"

    # Test Model indices
    @test isa(mod['A'], Chain)
    @test resname(mod['A']["H_20A"]) == "VAL"
    @test_throws KeyError mod['C']
    @test_throws MethodError mod["A"]

    # Test Model iteration
    model_list = collect(mod)
    @test isa(model_list, Vector{Chain})
    @test length(model_list) == 2
    @test chainid(model_list[2]) == 'B'


    # Test ProteinStructure getters/setters
    @test structurename(struc) == "Test structure"
    @test modelnumbers(struc) == [1, 3]
    @test isa(models(struc), Dict{Int, Model})
    @test length(models(struc)) == 2
    @test resids(models(struc)[1]['A']) == ["10", "H_20A"]

    @test isa(defaultmodel(struc), Model)
    @test modelnumber(defaultmodel(struc)) == 1
    @test chainids(struc) == ['A', 'B']

    # Test ProteinStructure indices
    @test isa(struc[1], Model)
    @test isa(struc[3], Model)
    @test isa(struc['A'], Chain)
    @test_throws KeyError struc[2]
    @test_throws KeyError struc['C']
    @test_throws MethodError struc["A"]
    @test ishetres(struc[1]['A']["H_20A"])

    # Test ProteinStructure iteration
    struc_list = collect(struc)
    @test isa(struc_list, Vector{Model})
    @test length(struc_list) == 2
    @test modelnumber(struc_list[2]) == 3


    # Test selector functions
    ch_a = Chain('A')
    ch_a["10"] = Residue("ALA", 10, ' ', false, ch_a)
    res_a = ch_a["10"]
    ch_a["H_11"] = Residue("MG", 11, ' ', true, ch_a)
    res_b = ch_a["H_11"]
    ch_a["H_100"] = Residue("HOH", 100, ' ', true, ch_a)
    res_c = ch_a["H_100"]
    ch_a["10"]["CA"] = Atom(100, " CA ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_a)
    atom_a = ch_a["10"]["CA"]
    ch_a["H_11"]["MG"] = Atom(110, "MG", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_b)
    atom_b = ch_a["H_11"]["MG"]

    @test stdatomselector(atom_a)
    @test !stdatomselector(atom_b)
    @test !hetatomselector(atom_a)
    @test hetatomselector(atom_b)
    @test atomnameselector(atom_a, Set(["CA", "N", "C"]))
    @test atomnameselector(atom_a, ["CA", "N", "C"])
    @test !atomnameselector(atom_b, Set(["CA", "N", "C"]))
    @test !atomnameselector(atom_a, ["CA", "N", "C"], spaces=true)
    @test calphaselector(atom_a)
    @test !calphaselector(atom_b)
    @test !calphaselector(Atom(100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_b))
    @test backboneselector(atom_a)
    @test !backboneselector(atom_b)
    @test !backboneselector(Atom(100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_b))
    @test heavyatomselector(atom_a)
    @test !heavyatomselector(atom_b)
    @test !heavyatomselector(Atom(100, "H1", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " H", "  ", res_a))
    @test resnameselector(atom_a, Set(["ALA"]))
    @test resnameselector(atom_a, ["ALA"])
    @test !resnameselector(atom_b, Set(["ALA"]))
    @test resnameselector(res_a, Set(["ALA"]))
    @test !resnameselector(res_b, Set(["ALA"]))
    @test !waterselector(res_a)
    @test waterselector(res_c)
    @test stdresselector(res_a)
    @test !stdresselector(res_b)
    @test !hetresselector(res_a)
    @test hetresselector(res_b)
    @test !disorderselector(atom_a)
    @test disorderselector(dis_at)
    @test !disorderselector(res_a)
    @test disorderselector(dis_res)
    @test hydrogenselector(Atom(100, "H", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " H", "  ", res_a))
    @test !hydrogenselector(Atom(100, "H", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_a))
    @test hydrogenselector(Atom(100, "H1", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "  ", "  ", res_a))
    @test hydrogenselector(Atom(100, "1H", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "  ", "  ", res_a))
    @test !hydrogenselector(Atom(100, "NH1", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "  ", "  ", res_a))


    # Further tests for structural element ordering
    # Order when looping over a DisorderedAtom is the atom serial
    dis_at_ord = DisorderedAtom(Dict(
        'A' => Atom(102, "CA", 'A', [1.0, 2.0, 3.0], 0.3, 10.0, "C", ""),
        'B' => Atom(101, "CA", 'B', [1.0, 2.0, 3.0], 0.4, 10.0, "C", ""),
        'C' => Atom(100, "CA", 'C', [1.0, 2.0, 3.0], 0.3, 10.0, "C", ""),
    ), 'B')
    @test altlocids(dis_at_ord) == ['C', 'B', 'A']

    # Order when sorting an atom list is the atom serial
    atom_list_ord = AbstractAtom[
        DisorderedAtom(Dict(
            'A' => Atom(100, "CA", 'A', [1.0, 2.0, 3.0], 0.4, 10.0, "C", ""),
            'B' => Atom(104, "CA", 'B', [1.0, 2.0, 3.0], 0.6, 10.0, "C", "")
        ), 'B'),
        Atom(102, "CB", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", ""),
        Atom(103, "CG", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    ]
    @test map(atomname, sort(atom_list_ord)) == ["CB", "CG", "CA"]
    sort!(atom_list_ord)
    @test map(atomname, atom_list_ord) == ["CB", "CG", "CA"]

    # Order when sorting a residue list is chain ID, then stdres/hetres,
    # then residue number, then insertion code
    residues_ord = AbstractResidue[
        Residue("ALA", 201, 'A', false, Chain('A')),
        Residue("ALA", 203, ' ', false, Chain('A')),
        Residue("ALA", 200, ' ', true, Chain('A')),
        Residue("ALA", 201, 'B', false, Chain('A')),
        Residue("ALA", 202, ' ', false, Chain('A')),
        Residue("ALA", 300, ' ', false, Chain('B')),
        Residue("ALA", 201, ' ', true, Chain('A')),
        Residue("ALA", 201, ' ', false, Chain('A')),
        Residue("ALA", 201, 'A', true, Chain('A')),
        Residue("ALA", 100, ' ', false, Chain('B')),
        Residue("ALA", 203, ' ', true, Chain('A')),
        Residue("ALA", 200, ' ', false, Chain('A')),
    ]
    # Test resid
    @test map(resid, residues_ord) == ["201A", "203", "H_200", "201B", "202", "300", "H_201", "201", "H_201A", "100", "H_203", "200"]
    @test map(res -> resid(res, full=true), residues_ord) == ["201A:A", "203:A", "H_200:A", "201B:A", "202:A", "300:B", "H_201:A", "201:A", "H_201A:A", "100:B", "H_203:A", "200:A"]
    @test map(res -> resid(res, full=true), sort(residues_ord)) == ["200:A", "201:A", "201A:A", "201B:A", "202:A", "203:A", "H_200:A", "H_201:A", "H_201A:A", "H_203:A", "100:B", "300:B"]
    sort!(residues_ord)
    @test map(res -> resid(res, full=true), residues_ord) == ["200:A", "201:A", "201A:A", "201B:A", "202:A", "203:A", "H_200:A", "H_201:A", "H_201A:A", "H_203:A", "100:B", "300:B"]

    # Order of listing residue names in a DisorderedResidue is default then alphabetical
    dis_res_ord = DisorderedResidue(Dict(
        "THR" => Residue("THR", 201, ' ', false),
        "ALA" => Residue("ALA", 201, ' ', false),
        "ILE" => Residue("ILE", 201, ' ', false),
        "SER" => Residue("SER", 201, ' ', false),
        "VAL" => Residue("VAL", 201, ' ', false)
    ), "SER")
    @test defaultresname(dis_res_ord) == "SER"
    @test resnames(dis_res_ord) == ["SER", "ALA", "ILE", "THR", "VAL"]

    # Order when sorting chain IDs is character ordering with the empty chain ID at the end
    model_ord = Model(1, Dict(
        'A' => Chain('A'),
        ' ' => Chain(' '),
        '1' => Chain('1'),
        'a' => Chain('a'),
        'X' => Chain('X'),
    ), ProteinStructure())
    @test chainids(model_ord) == ['1', 'A', 'X', 'a', ' ']
end


@testset "Parsing" begin
    # Test parsevalue
    line = "ATOM     40  CB  LEU A   5      22.088  45.547  29.675  1.00 22.23           C  "
    @test parsevalue(line, 7, 11, Int) == 40
    @test parsevalue(line, 31, 38, Float64) == 22.088
    @test parsevalue(line, 13, 16, String) == " CB "
    @test parsevalue(line, 22, 22, Char) == 'A'
    @test_throws ErrorException parsevalue(line, 1, 4, Int)
    @test_throws ErrorException parsevalue(line, 1, 4, Bool)
    @test_throws ErrorException parsevalue(line, 79, 100, Int)

    # Test parsestrict
    line =   "ATOM    591  C   GLY A  80      29.876  54.131  35.806  1.00 40.97           C  "
    line_a = "ATOM    591  C   GLY A  80              54.131  35.806  1.00 40.97           C  "
    line_b = "ATOM    591  C   GLY A  80 "
    @test parsestrict(line, 7, 11, Int, "could not read atom serial number", 10) == 591
    @test parsestrict(line, 13, 16, String, "could not read atom name", 20) == " C  "
    @test_throws PDBParseError parsestrict(line_a, 31, 38, Float64, "could not read x coordinate", 10)
    @test_throws PDBParseError parsestrict(line_b, 31, 38, Float64, "could not read x coordinate", 10)
    @test_throws PDBParseError parsestrict(line, 7, 11, Bool, "could not read atom serial number", 10)

    # Test parselenient
    line =   "ATOM     40  CB  LEU A   5      22.088  45.547  29.675  1.00 22.23           C  "
    line_a = "ATOM     40  CB  LEU A   5      22.088  45.547  29.675  1.00 22.23              "
    line_b = "ATOM     40  CB  LEU A   5      22.088  45.547  29.675  1.00 22.23  "
    line_c = "ATOM     40  CB  LEU A   5      22.088  45.547  29.675       22.23           C  "
    @test parselenient(line, 77, 78, String, "  ") == " C"
    @test parselenient(line_a, 77, 78, String, "  ") == "  "
    @test parselenient(line_b, 77, 78, String, "  ") == "  "
    @test parselenient(line_b, 77, 78, String, " C") == " C"
    @test parselenient(line_c, 55, 60, Float64, 1.0) == 1.0
    @test parselenient(line, 77, 78, Bool, " N") == " N"

    # Test AtomRecord constructor
    line_a = "ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  "
    line_b = "HETATM 3474  O  B XX A 334A      8.802  62.000   8.672  1.00 39.15           O1-"
    line_c = "ATOM    669  CA  ILE A  90      xxxxxx  33.110  31.221  1.00 25.76           C  "
    line_d = "ATOM    669  CA  ILE A  90      31.743   "
    line_e = "REMARK   1 REFERENCE 1                                                          "
    at = AtomRecord(line_a, 10)
    @test !at.het_atom
    @test at.serial == 669
    @test at.atom_name == " CA "
    @test at.alt_loc_id == ' '
    @test at.res_name == "ILE"
    @test at.chain_id == 'A'
    @test at.res_number == 90
    @test at.ins_code == ' '
    @test at.coords == [31.743, 33.110, 31.221]
    @test at.occupancy == 1.00
    @test at.temp_fac == 25.76
    @test at.element == " C"
    @test at.charge == "  "
    at = AtomRecord(line_b)
    @test at.het_atom
    @test at.serial == 3474
    @test at.atom_name == " O  "
    @test at.alt_loc_id == 'B'
    @test at.res_name == " XX"
    @test at.chain_id == 'A'
    @test at.res_number == 334
    @test at.ins_code == 'A'
    @test at.coords == [8.802, 62.0, 8.672]
    @test at.occupancy == 1.00
    @test at.temp_fac == 39.15
    @test at.element == " O"
    @test at.charge == "1-"
    @test_throws PDBParseError AtomRecord(line_c)
    @test_throws PDBParseError AtomRecord(line_d)

    # Test parsing 1AKE (multiple chains, disordered atoms)
    struc = read(pdbfilepath("1AKE.pdb"), PDB)
    @test structurename(struc) == "1AKE.pdb"
    @test countmodels(struc) == 1
    @test modelnumbers(struc) == [1]
    @test countchains(struc) == 2
    @test countchains(struc[1]) == 2
    @test chainids(struc) == ['A', 'B']
    @test chainids(struc[1]) == ['A', 'B']
    @test resname(struc['A'][10]) == "GLY"
    @test !isdisorderedres(struc['A'][10])
    @test serial(struc['A'][200]["NZ"]) == 1555
    @test !isdisorderedatom(struc['A'][200]["NZ"])
    @test isdisorderedatom(struc['A'][167]["CD"])
    @test altlocids(struc['A'][167]["CD"]) == ['A', 'B']
    @test x(struc['A'][167]["CD"]) == 24.502
    @test x(struc['A'][167]["CD"]['A']) == 24.502
    @test x(struc['A'][167]["CD"]['B']) == 24.69

    """
    # Test collectatoms
    ats = collectatoms(struc)
    @test length(ats) == 3804
    @test isa(ats, Vector{AbstractAtom})
    @test isa(ats[70], Atom)
    @test isa(ats[1290], DisorderedAtom)
    @test serial(ats[1660]) == 3323
    ats = collectatoms(struc, hetatomselector)
    @test length(ats) == 492
    @test serial(ats[60]) == 3443
    ats = collectatoms(struc, disorderselector)
    @test length(ats) == 12
    @test serial(ats[10]) == 3338
    ats = collectatoms(struc, stdatomselector, disorderselector)
    @test length(ats) == 5
    @test serial(ats[4]) == 1294
    ats = collectatoms(struc[1])
    @test length(ats) == 3804
    @test serial(ats[1660]) == 3323
    ats = collectatoms(struc['A'])
    @test length(ats) == 1954
    @test serial(ats[240]) == 240
    ats = collectatoms(struc['A'][50])
    @test length(ats) == 9
    @test serial(ats[4]) == 358
    ats = collectatoms(struc['A'][50]["CA"])
    @test length(ats) == 1
    @test isa(ats, Vector{AbstractAtom})
    @test serial(ats[1]) == 356
    ats = collectatoms(struc['A'][167]["CZ"])
    @test length(ats) == 1
    @test isa(ats, Vector{AbstractAtom})
    @test isa(ats[1], DisorderedAtom)
    @test serial(ats[1]) == 1292
    ats = collectatoms(Chain[struc['B'], struc['A']])
    @test length(ats) == 3804
    @test serial(ats[5]) == 1667
    ats = collectatoms(Residue[struc['A'][51], struc['A'][50]])
    @test length(ats) == 17
    @test serial(ats[10]) == 356
    ats = collectatoms(AbstractResidue[struc['A'][51], struc['A'][50]])
    @test length(ats) == 17
    @test serial(ats[10]) == 356
    ats = collectatoms(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test length(ats) == 2
    @test isa(ats, Vector{Atom})
    @test serial(ats[1]) == 356
    ats = collectatoms(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(ats) == 2
    @test isa(ats, Vector{DisorderedAtom})
    @test serial(ats[2]) == 1292
    ats = collectatoms(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(ats) == 2
    @test isa(ats, Vector{AbstractAtom})
    @test serial(ats[2]) == 1292


    # Test countatoms
    @test countatoms(struc) == 3804
    @test countatoms(struc[1]) == 3804
    @test countatoms(struc['A']) == 1954
    @test countatoms(struc['A'][50]) == 9
    @test countatoms(struc['A'][50]["CA"]) == 1
    @test countatoms([struc['A'], struc['B']]) == 3804
    @test countatoms(AbstractResidue[struc['A'][50], struc['A'][51]]) == 17
    @test countatoms(Residue[struc['A'][50], struc['A'][51]]) == 17
    @test countatoms(collectatoms(struc['A'])) == 1954
    @test countatoms(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]]) == 2
    @test countatoms(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]]) == 2

    @test countatoms(struc['A'], stdatomselector) == 1656
    @test countatoms(struc['A'], hetatomselector) == 298
    @test countatoms(struc['A'], stdatomselector, disorderselector) == 5

    @test countatoms(ProteinStructure()) == 0
    @test countatoms(Model()) == 0
    @test countatoms(Chain('X')) == 0
    @test countatoms(Residue("ALA", 'A', 100, ' ', false)) == 0


    # Test collectresidues
    residues = collectresidues(struc)
    @test length(residues) == 808
    @test isa(residues, Vector{AbstractResidue})
    @test isa(residues[50], Residue)
    @test resnumber(residues[220]) == 305
    residues = collectresidues(struc, hetresselector)
    @test length(residues) == 380
    @test resnumber(residues[370]) == 725
    residues = collectresidues(struc, stdresselector, res -> chainid(res) == 'A')
    @test length(residues) == 214
    @test resnumber(residues[200]) == 200
    residues = collectresidues(struc[1])
    @test length(residues) == 808
    @test resnumber(residues[220]) == 305
    residues = collectresidues(struc['A'])
    @test length(residues) == 456
    @test resnumber(residues[220]) == 305
    residues = collectresidues(struc['A'][50])
    @test length(residues) == 1
    @test isa(residues, Vector{AbstractResidue})
    @test resnumber(residues[1]) == 50
    residues = collectresidues(struc['A'][50]["CA"])
    @test length(residues) == 1
    @test isa(residues, Vector{AbstractResidue})
    @test isa(residues[1], Residue)
    @test resnumber(residues[1]) == 50
    residues = collectresidues(struc['A'][167]["CZ"])
    @test length(residues) == 1
    @test isa(residues, Vector{AbstractResidue})
    @test isa(residues[1], Residue)
    @test resnumber(residues[1]) == 167
    residues = collectresidues(Chain[struc['B'], struc['A']])
    @test length(residues) == 808
    @test resid(residues[5]; full=true) == "5:B"
    residues = collectresidues(Residue[struc['A'][51], struc['A'][50]])
    @test length(residues) == 2
    @test isa(residues, Vector{Residue})
    @test resnumber(residues[1]) == 50
    residues = collectresidues(AbstractResidue[struc['A'][51], struc['A'][50]])
    @test length(residues) == 2
    @test isa(residues, Vector{AbstractResidue})
    @test resnumber(residues[1]) == 50
    residues = collectresidues(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test length(residues) == 2
    @test isa(residues, Vector{AbstractResidue})
    @test resnumber(residues[2]) == 51
    residues = collectresidues(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(residues) == 1
    @test isa(residues, Vector{AbstractResidue})
    @test atomnames(residues[1]) == ["CD", "CZ"]
    residues = collectresidues(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(residues) == 2
    @test isa(residues, Vector{AbstractResidue})
    @test atomnames(residues[1]) == ["CA"]


    # Test countresidues
    @test countresidues(struc) == 808
    @test countresidues(struc[1]) == 808
    @test countresidues(struc['A']) == 456
    @test countresidues(struc['A'][50]) == 1
    @test countresidues(struc['A'][50]["CA"]) == 1
    @test countresidues([struc['A'], struc['B']]) == 808
    @test countresidues(AbstractResidue[struc['A'][50], struc['A'][51]]) == 2
    @test countresidues(Residue[struc['A'][50], struc['A'][51]]) == 2
    @test countresidues(collectatoms(struc['A'])) == 456
    @test countresidues(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]]) == 1
    @test countresidues(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]]) == 2

    @test countresidues(struc['A'], stdresselector) == 214
    @test countresidues(struc['A'], hetresselector) == 242
    @test countresidues(struc, stdresselector, res -> chainid(res) == 'A') == 214

    @test countresidues(ProteinStructure()) == 0
    @test countresidues(Model()) == 0
    @test countresidues(Chain('X')) == 0
    @test countresidues(Residue("ALA", 'A', 100, ' ', false)) == 1


    # Test formatomlist
    atom_a = Atom(false, 103, "N", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "N", "")
    atom_b = Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.4, 10.0, "C", "")
    atom_c = Atom(false, 102, "C", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    atom_d = Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.6, 10.0, "C", "")
    atom_list = formatomlist([atom_a, atom_b, atom_c, atom_d])
    @test length(atom_list) == 3
    @test isa(atom_list[1], DisorderedAtom)
    @test isa(atom_list[2], Atom)
    @test isa(atom_list[3], Atom)
    @test map(length, atom_list) == [2, 1, 1]
    @test map(serial, atom_list) == [101, 102, 103]
    @test defaultaltlocid(atom_list[1]) == 'B'
    atom_list = formatomlist([atom_a, atom_b, atom_c, atom_d]; remove_disorder=true)
    @test length(atom_list) == 3
    @test map(x -> isa(x, Atom), atom_list) == [true, true, true]
    @test map(length, atom_list) == [1, 1, 1]
    @test map(serial, atom_list) == [101, 102, 103]
    atom_a = Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    atom_b = Atom(false, 101, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    @test_throws ErrorException formatomlist([atom_a, atom_b])

    # Test by unfolding 1AKE atoms and re-forming them
    atoms_unfold = Atom[]
    for atom in collectatoms(struc['A'])
        append!(atoms_unfold, collect(atom))
    end
    atom_list = formatomlist(atoms_unfold)
    @test length(atom_list) == 1954
    @test sum(map(x -> isa(x, Atom), atom_list)) == 1942
    @test isa(atom_list[1289], DisorderedAtom)
    @test length(atom_list[1289]) == 2
    @test serial(atom_list[1934]) == 3661
    atom_list = formatomlist(atoms_unfold; remove_disorder=true)
    @test length(atom_list) == 1954
    @test sum(map(x -> isa(x, Atom), atom_list)) == 1954
    @test serial(atom_list[1934]) == 3661


    # Test choosedefaultaltlocid
    atom_a = Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.4, 10.0, "C", "")
    atom_b = Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.6, 10.0, "C", "")
    @test choosedefaultaltlocid(atom_a, atom_b) == 'B'
    @test choosedefaultaltlocid(atom_b, atom_a) == 'B'
    atom_a = Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 10.0, "C", "")
    atom_b = Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 10.0, "C", "")
    @test choosedefaultaltlocid(atom_a, atom_b) == 'A'
    @test choosedefaultaltlocid(atom_b, atom_a) == 'A'


    # Test applyselectors
    atoms = collectatoms(struc)
    # Not providing any selector functions just returns the input list
    atoms_min = applyselectors(atoms)
    @test length(atoms_min) == length(atoms)
    @test map(serial, atoms_min) == map(serial, atoms)
    applyselectors!(atoms_min)
    @test length(atoms_min) == length(atoms)
    @test map(serial, atoms_min) == map(serial, atoms)
    atoms_min = applyselectors(atoms, stdatomselector)
    @test length(atoms_min) == 3312
    @test serial(atoms_min[2000]) == 2006
    applyselectors!(atoms, stdatomselector)
    @test length(atoms) == 3312
    @test serial(atoms[2000]) == 2006
    atoms = collectatoms(struc)
    atoms_min = applyselectors(atoms, stdatomselector, disorderselector)
    @test length(atoms_min) == 5
    @test serial(atoms_min[4]) == 1294
    applyselectors!(atoms, stdatomselector, disorderselector)
    @test length(atoms) == 5
    @test serial(atoms[4]) == 1294

    residues = collectresidues(struc)
    residues_min = applyselectors(residues)
    @test length(residues_min) == length(residues)
    @test map(res -> resid(res; full=true), residues_min) == map(res -> resid(res; full=true), residues)
    applyselectors!(residues_min)
    @test length(residues_min) == length(residues)
    @test map(res -> resid(res; full=true), residues_min) == map(res -> resid(res; full=true), residues)
    residues_min = applyselectors(residues, waterselector)
    @test length(residues_min) == 378
    @test resid(residues_min[300]; full=true) == "H_657:B"
    applyselectors!(residues, waterselector)
    @test length(residues) == 378
    @test resid(residues[300]; full=true) == "H_657:B"
    residues = collectresidues(struc)
    # Test anonymous selector function
    residues_min = applyselectors(residues, stdresselector, res -> chainid(res) == 'A')
    @test length(residues_min) == 214
    @test resid(residues_min[200]; full=true) == "200:A"
    applyselectors!(residues, stdresselector, res -> chainid(res) == 'A')
    @test length(residues) == 214
    @test resid(residues[200]; full=true) == "200:A"


    # Test parsing options
    struc = read(pdbfilepath("1AKE.pdb"), PDB; structure_name="New name")
    @test structurename(struc) == "New name"
    @test countatoms(struc) == 3804

    struc = read(pdbfilepath("1AKE.pdb"), PDB; read_het_atoms=false)
    @test countatoms(struc) == 3312
    @test serial(collectatoms(struc)[2000]) == 2006
    @test sum(map(ishetatom, collectatoms(struc))) == 0

    struc = read(pdbfilepath("1AKE.pdb"), PDB; read_std_atoms=false)
    @test countatoms(struc) == 492
    @test serial(collectatoms(struc)[400]) == 3726
    @test sum(map(ishetatom, collectatoms(struc))) == 492

    struc = read(pdbfilepath("1AKE.pdb"), PDB; read_het_atoms=false, read_std_atoms=false)
    @test countatoms(struc) == 0
    @test countresidues(struc) == 0

    struc = read(pdbfilepath("1AKE.pdb"), PDB; remove_disorder=true)
    @test countatoms(struc) == 3804
    @test sum(map(isdisorderedatom, collectatoms(struc))) == 0
    @test tempfac(struc['A'][167]["NE"]) == 23.32

    struc = read(pdbfilepath("1AKE.pdb"), PDB, backboneselector)
    @test countatoms(struc) == 1284
    @test countatoms(struc, backboneselector) == 1284
    @test serial(collectatoms(struc)[1000]) == 2566

    struc = read(pdbfilepath("1AKE.pdb"), PDB, stdatomselector, disorderselector)
    @test countatoms(struc) == 5
    @test sum(map(isdisorderedatom, collectatoms(struc))) == 5
    @test sum(map(ishetatom, collectatoms(struc))) == 0
    @test serial(struc['A'][167]["CZ"]) == 1292


    # Test parsing from stream
    open(pdbfilepath("1AKE.pdb"), "r") do file
        struc = read(file, PDB)
        @test countatoms(struc) == 3804
        @test countresidues(struc) == 808
    end

    # Test parsing 1EN2 (disordered residue)
    struc = read(pdbfilepath("1EN2.pdb"), PDB)
    @test modelnumbers(struc) == [1]
    @test chainids(struc[1]) == ['A']
    @test serial(struc['A'][48]["CA"]) == 394
    @test isdisorderedres(struc['A'][10])
    # Note default is first added, not highest occupancy
    @test defaultresname(struc['A'][10]) == "SER"
    @test resname(disorderedres(struc['A'][10], "SER")) == "SER"
    @test resname(disorderedres(struc['A'][10], "GLY")) == "GLY"
    # Atoms in a disordered residue are not necessarily disordered
    @test !isdisorderedatom(struc['A'][10]["CA"])
    @test altlocid(struc['A'][10]["CA"]) == 'A'
    @test isdisorderedres(struc['A'][16])
    @test defaultresname(struc['A'][16]) == "ARG"
    @test resname(disorderedres(struc['A'][16], "ARG")) == "ARG"
    @test resname(disorderedres(struc['A'][16], "TRP")) == "TRP"
    @test !isdisorderedatom(disorderedres(struc['A'][16], "TRP")["CA"])
    @test altlocid(disorderedres(struc['A'][16], "TRP")["CA"]) == 'C'
    @test isdisorderedatom(struc['A'][16]["CA"])
    @test altlocids(struc['A'][16]["CA"]) == ['A', 'B']
    @test defaultaltlocid(struc['A'][16]["CA"]) == 'A'
    @test occupancy(struc['A'][16]["CA"]) == 0.22
    atoms = collectatoms(DisorderedResidue[struc['A'][10], struc['A'][16]])
    @test length(atoms) == 17
    @test isa(atoms, Vector{AbstractAtom})
    @test serial(atoms[10]) == 113
    @test isa(atoms[10], DisorderedAtom)
    @test countatoms(DisorderedResidue[struc['A'][10], struc['A'][16]]) == 17
    residues = collectresidues(DisorderedResidue[struc['A'][16], struc['A'][10]])
    @test length(residues) == 2
    @test isa(residues, Vector{DisorderedResidue})
    @test resnumber(residues[1]) == 10
    @test countresidues(DisorderedResidue[struc['A'][16], struc['A'][10]]) == 2


    # Test parsing 1SSU (multiple models)
    struc = read(pdbfilepath("1SSU.pdb"), PDB)
    # Test countmodels
    @test countmodels(struc) == 20
    @test modelnumbers(struc) == collect(1:20)
    @test countchains(struc) == 1
    @test countchains(struc[5]) == 1
    @test chainids(struc) == ['A']
    @test chainids(struc[5]) == ['A']
    @test serial(struc[10]['A'][40]["HB2"]) == 574
    # Note counting and collecting on model 1 by default
    @test countatoms(struc) == 756
    @test map(countatoms, [struc[i] for i in modelnumbers(struc)]) == 756 * ones(Int, 20)
    @test countatoms(struc, hydrogenselector) == 357
    atoms = collectatoms(Model[struc[5], struc[10]])
    @test length(atoms) == 1512
    @test z(atoms[20]) == -14.782
    @test z(atoms[1000]) == -3.367
    @test countatoms(Model[struc[5], struc[10]]) == 1512
    atoms_raw = Atom[atom for atom in atoms]
    @test_throws ErrorException formatomlist(atoms_raw)
    residues = collectresidues(Model[struc[5], struc[10]])
    @test length(residues) == 102
    @test y(residues[10]["O"]) == -1.612
    @test y(residues[100]["O"]) == -13.184
    @test countresidues(Model[struc[5], struc[10]]) == 102


    # Test organise
    struc_new = organise(Model[struc[5], struc[3]])
    @test isa(struc_new, ProteinStructure)
    @test structurename(struc_new) == ""
    @test modelnumbers(struc_new) == [3, 5]
    @test countatoms(struc_new[3]) == 756
    @test countatoms(struc_new[5]) == 756
    @test_throws KeyError struc_new[1]
    struc_new = organise(Model[struc[5], struc[3]]; structure_name="new struc")
    @test structurename(struc_new) == "new struc"
    struc_new = organise(struc[5])
    @test isa(struc_new, ProteinStructure)
    @test modelnumbers(struc_new) == [5]
    @test countatoms(struc_new[5]) == 756

    struc = read(pdbfilepath("1AKE.pdb"), PDB)
    model_new = organise(Chain[struc['B'], struc['A']])
    @test isa(model_new, Model)
    @test modelnumber(model_new) == 1
    @test chainids(model_new) == ['A', 'B']
    @test countatoms(model_new['A']) == 1954
    @test countatoms(model_new['B']) == 1850
    @test_throws KeyError model_new[' ']
    model_new = organise(Chain[struc['B'], struc['A']]; model_number=100)
    @test modelnumber(model_new) == 100
    model_new = organise(struc['B'])
    @test isa(model_new, Model)
    @test chainids(model_new) == ['B']
    @test countatoms(model_new['B']) == 1850

    chains_new = organise(shuffle(collectresidues(struc)))
    @test isa(chains_new, Vector{Chain})
    @test length(chains_new) == 2
    @test chainid(chains_new[1]) == 'A'
    @test map(chainid, chains_new[2]) == ['B' for i in 1:352]
    @test x(chains_new[2]["H_725"]["O"]) == 34.939
    chains_new = organise(Residue[struc['A'][10], struc['A'][11], struc['B'][10]])
    @test isa(chains_new, Vector{Chain})
    @test length(chains_new) == 2
    @test chainid(chains_new[2]) == 'B'
    @test x(chains_new[2][10]["C"]) == 23.612
    struc = read(pdbfilepath("1EN2.pdb"), PDB)
    chains_new = organise(DisorderedResidue[struc['A'][10], struc['A'][16]])
    @test isa(chains_new, Vector{Chain})
    @test length(chains_new) == 1
    @test chainid(chains_new[1]) == 'A'
    @test x(chains_new[1][10]["C"]) == -5.157
    struc = read(pdbfilepath("1AKE.pdb"), PDB)
    chains_new = organise(struc['A'][50])
    @test isa(chains_new, Vector{Chain})
    @test length(chains_new) == 1
    @test chainid(chains_new[1]) == 'A'
    @test x(chains_new[1]["50"]["NZ"]) == 36.415

    residues_new = organise(shuffle(collectatoms(struc)))
    @test isa(residues_new, Vector{AbstractResidue})
    @test length(residues_new) == 808
    @test map(chainid, residues_new) == [['A' for i in 1:456]; ['B' for i in 1:352]]
    @test tempfac(residues_new[10]["C"]) == 19.36
    @test tempfac(residues_new[215]["PB"]) == 16.65
    residues_new = organise(Atom[struc['A'][50]["NZ"], struc['A'][51]["OD2"]])
    @test isa(residues_new, Vector{AbstractResidue})
    @test length(residues_new) == 2
    @test resid(residues_new[2]) == "51"
    @test countatoms(residues_new[1]) == 1
    residues_new = organise(DisorderedAtom[struc['A'][167]["NH1"], struc['A'][167]["NH2"]])
    @test isa(residues_new, Vector{AbstractResidue})
    @test length(residues_new) == 1
    @test resid(residues_new[1]) == "167"
    @test countatoms(residues_new[1]) == 2
    residues_new = organise(struc['A'][50]["NZ"])
    @test isa(residues_new, Vector{AbstractResidue})
    @test length(residues_new) == 1
    @test resid(residues_new[1]) == "50"
    @test countatoms(residues_new[1]) == 1

    struc = read(pdbfilepath("1EN2.pdb"), PDB)
    atoms = [
        collect(struc['A'][9])
        collect(disorderedres(struc['A'][10], "SER"))
        collect(disorderedres(struc['A'][10], "GLY"))
    ]
    residues_new = organise(atoms)
    @test isa(residues_new, Vector{AbstractResidue})
    @test length(residues_new) == 2
    @test isa(residues_new[2], DisorderedResidue)
    @test isa(residues_new[1], Residue)
    @test defaultresname(residues_new[2]) == "SER"
    @test countatoms(disorderedres(residues_new[2], "SER")) == 6
    @test countatoms(disorderedres(residues_new[2], "GLY")) == 4
    @test !isdisorderedatom(residues_new[2]["CB"])


    # Test organisemodel
    struc = read(pdbfilepath("1AKE.pdb"), PDB)
    model_new = organisemodel(Chain[struc['B'], struc['A']])
    @test isa(model_new, Model)
    @test modelnumber(model_new) == 1
    @test chainids(model_new) == ['A', 'B']
    @test countatoms(model_new['A']) == 1954
    model_new = organisemodel(Chain[struc['B'], struc['A']]; model_number=5)
    @test modelnumber(model_new) == 5
    model_new = organisemodel(struc['B'])
    @test isa(model_new, Model)
    @test chainids(model_new) == ['B']
    @test countatoms(model_new['B']) == 1850

    model_new = organisemodel(shuffle(collectresidues(struc)))
    @test isa(model_new, Model)
    @test chainids(model_new) == ['A', 'B']
    @test countatoms(model_new['A']) == 1954
    model_new = organisemodel(struc['A'][50])
    @test isa(model_new, Model)
    @test chainids(model_new) == ['A']
    @test countatoms(model_new) == 9

    model_new = organisemodel(shuffle(collectatoms(struc)))
    @test isa(model_new, Model)
    @test chainids(model_new) == ['A', 'B']
    @test countatoms(model_new['A']) == 1954
    model_new = organisemodel(struc['A'][50]["NZ"])
    @test isa(model_new, Model)
    @test chainids(model_new) == ['A']
    @test countatoms(model_new) == 1


    # Test organisestructure
    struc = read(pdbfilepath("1SSU.pdb"), PDB)
    struc_new = organisestructure(Model[struc[5], struc[3]])
    @test isa(struc_new, ProteinStructure)
    @test structurename(struc_new) == ""
    @test modelnumbers(struc_new) == [3, 5]
    @test countatoms(struc_new[3]) == 756
    struc_new = organisestructure(Model[struc[5], struc[3]]; structure_name="new struc")
    @test structurename(struc_new) == "new struc"
    struc_new = organisestructure(struc[5])
    @test isa(struc_new, ProteinStructure)
    @test modelnumbers(struc_new) == [5]
    @test countatoms(struc_new[5]) == 756

    struc = read(pdbfilepath("1AKE.pdb"), PDB)
    struc_new = organisestructure(Chain[struc['B'], struc['A']])
    @test isa(struc_new, ProteinStructure)
    @test modelnumbers(struc_new) == [1]
    @test modelnumber(defaultmodel(struc_new)) == 1
    @test countatoms(struc_new) == 3804
    struc_new = organisestructure(Chain[struc['B'], struc['A']]; model_number=7)
    @test modelnumbers(struc_new) == [7]
    @test modelnumber(defaultmodel(struc_new)) == 7
    struc_new = organisestructure(struc['B'])
    @test isa(struc_new, ProteinStructure)
    @test modelnumbers(struc_new) == [1]
    @test countatoms(struc_new) == 1850

    struc_new = organisestructure(shuffle(collectresidues(struc)))
    @test isa(struc_new, ProteinStructure)
    @test modelnumbers(struc_new) == [1]
    @test countatoms(struc_new) == 3804
    struc_new = organisestructure(struc['A'][50])
    @test isa(struc_new, ProteinStructure)
    @test modelnumbers(struc_new) == [1]
    @test countatoms(struc_new) == 9

    struc_new = organisestructure(shuffle(collectatoms(struc)))
    @test isa(struc_new, ProteinStructure)
    @test modelnumbers(struc_new) == [1]
    @test countatoms(struc_new) == 3804
    struc_new = organisestructure(struc['A'][50]["NZ"])
    @test isa(struc_new, ProteinStructure)
    @test modelnumbers(struc_new) == [1]
    @test countatoms(struc_new) == 1
    """


    # Test parser error handling
    error = PDBParseError("message", 10, "line")
    showerror(DevNull, error)
    # Missing coordinate (blank string)
    @test_throws PDBParseError read(pdbfilepath("1AKE_err_a.pdb"), PDB)
    # Missing chain ID (line ends early)
    @test_throws PDBParseError read(pdbfilepath("1AKE_err_b.pdb"), PDB)
    # Bad MODEL record
    @test_throws PDBParseError read(pdbfilepath("1SSU_err.pdb"), PDB)
    # Duplicate atom names in same residue
    @test_throws ErrorException read(pdbfilepath("1AKE_err_c.pdb"), PDB)
    # Non-existent file
    @test_throws SystemError read(pdbfilepath("non_existent_file.pdb"), PDB)
end


@testset "Writing" begin
    # Test spacestring
    @test spacestring(1.5, 5) == "  1.5"
    @test spacestring("A", 3) == "  A"
    @test spacestring('A', 3) == "  A"
    @test_throws AssertionError spacestring(1.456789, 5)
    @test_throws AssertionError spacestring("ABCDEF", 3)


    # Test spaceatomname
    @test spaceatomname(Atom(1, " CA ", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " C", "  ")) == " CA "
    @test spaceatomname(Atom(1, "N",    ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " N", "  ")) == " N  "
    @test spaceatomname(Atom(1, "N",    ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "  ", "  ")) == " N  "
    @test spaceatomname(Atom(1, "CA",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " C", "  ")) == " CA "
    @test spaceatomname(Atom(1, "NE1",  ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " N", "  ")) == " NE1"
    @test spaceatomname(Atom(1, "2HD1", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ")) == "2HD1"
    @test spaceatomname(Atom(1, "HH11", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ")) == "HH11"
    @test spaceatomname(Atom(1, "1H",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ")) == "1H  "
    @test spaceatomname(Atom(1, "MG",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "MG", "  ")) == "MG  "
    @test spaceatomname(Atom(1, "MG",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "  ", "  ")) == " MG "
    @test_throws AssertionError spaceatomname(Atom(1, "11H",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  "))
    @test_throws AssertionError spaceatomname(Atom(1, "11H11", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  "))
    @test_throws AssertionError spaceatomname(Atom(1, "1MG",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "MG", "  "))


    # Test pdbline
    ch_a = Chain('A')
    ch_a["1"] = Residue("ALA", 1, ' ', false, ch_a)
    ch_a["1"][" N  "] = Atom(10, " N  ", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " N", "  ", ch_a["1"])
    line = join(pdbline(ch_a["1"][" N  "]))
    @test line == "ATOM     10  N   ALA A   1         0.0     0.0     0.0   1.0   0.0           N  "
    ch_b = Chain('B')
    ch_b["H_20"] = Residue("X", 20, ' ', true, ch_b)
    ch_b["H_20"]["C"] = Atom(101, "C", 'A', [10.5, 20.12345, -5.1227], 0.50, 50.126, "C", "1+", ch_b["H_20"])
    line = join(pdbline(ch_b["H_20"]["C"]))
    @test line == "HETATM  101  C  A  X B  20        10.5  20.123  -5.123   0.5 50.13           C1+"
    ch_b["H_20"]["11H11"] = Atom(1, "11H11", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ")
    @test_throws AssertionError pdbline(ch_b["H_20"]["11H11"])


    # Test writepdb and writepdblines
    # Counts lines in a file
    function countlines(filename::AbstractString)
        counter = 0
        open(filename, "r") do file
            for line in eachline(file)
                counter += 1
            end
        end
        return counter
    end

    struc = read(pdbfilepath("1SSU.pdb"), PDB)
    # All writing is done to one temporary file which is removed at the end
    temp_filename = tempname()
    writepdb(temp_filename, struc)
    @test countlines(temp_filename) == 15160
    struc_written = read(temp_filename, PDB)
    @test isa(struc_written, ProteinStructure)
    @test modelnumbers(struc_written) == collect(1:20)
    @test countatoms(struc_written) == 756
    @test z(struc_written[4]['A']["30"]["OG"]) == -2.177
    @test atomnames(struc_written[15]['A']["39"]) == ["N", "CA", "C", "O", "CB", "SG", "H", "HA", "HB2", "HB3"]

    # Test writing to stream
    open(temp_filename, "w") do file
        writepdb(file, struc)
    end
    @test countlines(temp_filename) == 15160
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == collect(1:20)
    @test countatoms(struc_written) == 756
    @test z(struc_written[4]['A']["30"]["OG"]) == -2.177
    @test atomnames(struc_written[15]['A']["39"]) == ["N", "CA", "C", "O", "CB", "SG", "H", "HA", "HB2", "HB3"]

    open(temp_filename, "w") do file
        writepdblines(file, struc)
    end
    @test countlines(temp_filename) == 756
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 756
    @test resname(struc_written['A'][13]["CE1"]) == "PHE"

    # Test selectors
    struc = read(pdbfilepath("1AKE.pdb"), PDB)
    writepdb(temp_filename, struc, hetatomselector)
    @test countlines(temp_filename) == 499
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 492
    @test chainids(struc_written) == ['A', 'B']
    @test tempfac(struc_written['B']["H_705"]["O"]) == 64.17
    writepdb(temp_filename, struc, stdatomselector, disorderselector)
    @test countlines(temp_filename) == 10
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 5
    @test sum(map(isdisorderedatom, collectatoms(struc_written))) == 5
    @test defaultaltlocid(struc_written['A'][167]["NH1"]) == 'A'


    """
    # Test writing different element types
    writepdb(temp_filename, struc[1])
    @test countlines(temp_filename) == 3816
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 3804
    writepdb(temp_filename, struc['A'])
    @test countlines(temp_filename) == 1966
    struc_written = read(temp_filename, PDB)
    @test chainids(struc_written) == ['A']
    writepdb(temp_filename, struc['A'][50])
    @test countlines(temp_filename) == 9
    struc_written = read(temp_filename, PDB)
    @test chainids(struc_written) == ['A']
    @test countresidues(struc_written) == 1
    @test atomnames(struc_written['A'][50]) == ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"]
    writepdb(temp_filename, struc['A'][50]["CA"])
    @test countlines(temp_filename) == 1
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 1
    @test !isdisorderedatom(collectatoms(struc_written)[1])
    @test serial(collectatoms(struc_written)[1]) == 356
    writepdb(temp_filename, struc['A'][167]["CZ"])
    @test countlines(temp_filename) == 2
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 1
    @test isdisorderedatom(collectatoms(struc_written)[1])
    @test length(collectatoms(struc_written)[1]) == 2
    @test tempfac(collectatoms(struc_written)[1]) == 16.77
    writepdb(temp_filename, Chain[struc['A'], struc['B']])
    @test countlines(temp_filename) == 3816
    struc_written = read(temp_filename, PDB)
    @test chainids(struc_written) == ['A', 'B']
    @test countatoms(struc_written['A']) == 1954
    @test countatoms(struc_written['B']) == 1850
    @test altlocids(struc_written['A']["H_215"]["O1G"]) == ['A', 'B']
    writepdb(temp_filename, AbstractResidue[struc['A'][51], struc['A'][50]])
    @test countlines(temp_filename) == 17
    struc_written = read(temp_filename, PDB)
    @test countresidues(struc_written) == 2
    @test map(resid, collectresidues(struc_written)) == ["50", "51"]
    @test countatoms(struc_written) == 17
    writepdb(temp_filename, AbstractAtom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test countlines(temp_filename) == 2
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 2
    @test !ishetatom(struc_written['A'][51]["CA"])

    # Test multiple model writing
    struc = read(pdbfilepath("1SSU.pdb"), PDB)
    writepdb(temp_filename, Model[struc[10], struc[5]])
    @test countlines(temp_filename) == 1516
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == [5, 10]
    @test modelnumber(defaultmodel(struc_written)) == 5
    @test countatoms(struc_written[5]) == 756
    @test countatoms(struc_written[10]) == 756
    @test_throws KeyError struc_written[1]

    # Test disordered residue writing
    struc = read(pdbfilepath("1EN2.pdb"), PDB)
    writepdb(temp_filename, struc)
    @test countlines(temp_filename) == 819
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 754
    @test isa(struc_written['A'][15], Residue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][16]) == "ARG"
    @test isa(struc_written['A'][16]["N"], DisorderedAtom)
    @test defaultaltlocid(struc_written['A'][16]["N"]) == 'A'
    @test isa(disorderedres(struc_written['A'][16], "TRP")["N"], Atom)
    @test countatoms(struc_written['A'][16]) == 11
    @test countatoms(disorderedres(struc_written['A'][16], "TRP")) == 14
    writepdb(temp_filename, AbstractResidue[struc['A'][16], struc['A'][10]])
    @test countlines(temp_filename) == 46
    struc_written = read(temp_filename, PDB)
    @test countresidues(struc_written) == 2
    @test isa(struc_written['A'][10], DisorderedResidue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][10]) == "SER"
    @test isa(disorderedres(struc_written['A'][10], "GLY")["O"], Atom)
    @test altlocid(disorderedres(struc_written['A'][10], "GLY")["O"]) == 'B'
    @test countatoms(struc_written['A'][10]) == 6
    @test countatoms(struc_written['A'][16]) == 11

    @test_throws AssertionError writepdb(temp_filename, Atom(false, 1, "11H11", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "H", ""))
    """

    rm(temp_filename)
end


@testset "Spatial" begin
    # Test coordarray
    atom = Atom(100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ")
    cs = coordarray(atom)
    @test size(cs) == (3,1)
    @test cs[1] == 1.0
    @test cs[2] == 2.0
    @test cs[3] == 3.0

    struc_1AKE = read(pdbfilepath("1AKE.pdb"), PDB)
    cs = coordarray(struc_1AKE)
    @test size(cs) == (3,3804)
    @test cs[1,3787] == 20.135
    @test cs[2,3787] == -10.789
    @test cs[3,3787] == -1.732
    cs = coordarray(struc_1AKE['A'], calphaselector)
    @test size(cs) == (3,214)
    @test cs[1,10] == 17.487
    @test cs[2,10] == 42.426
    @test cs[3,10] == 19.756
    @test coordarray(cs) == cs


    # Test rmsd
    cs_one = [
        0.0 0.0
        1.0 0.0
        1.0 3.0
    ]
    cs_two = [
        0.0 2.0
        2.0 0.0
        1.0 3.0
    ]
    # This line gives a @simd warning when --inline=no
    @test isapprox(rmsd(cs_one, cs_two), sqrt(5/2))
    cs_one = [
        0.0 0.0 1.0
        1.0 0.0 2.0
        1.0 3.0 3.0
    ]
    cs_two = [
        0.0 2.0
        2.0 0.0
        1.0 3.0
    ]
    @test_throws AssertionError rmsd(cs_one, cs_two)

    struc_1SSU = read(pdbfilepath("1SSU.pdb"), PDB)
    @test isapprox(rmsd(struc_1SSU[1], struc_1SSU[2], calphaselector), 4.1821925809691889)
    @test isapprox(rmsd(struc_1SSU[5], struc_1SSU[6], backboneselector), 5.2878196391279939)
    @test_throws AssertionError rmsd(struc_1SSU[1]['A'][8], struc_1SSU[1]['A'][9])


    # Test displacements
    cs_one = [
        0.0 0.0
        1.0 0.0
        1.0 3.0
    ]
    cs_two = [
        0.0 2.0
        2.0 0.0
        1.0 4.0
    ]
    # This line gives a @simd warning when --inline=no
    @test isapprox(displacements(cs_one, cs_two), [1.0, sqrt(5)])
    cs_one = [
        0.0 0.0 1.0
        1.0 0.0 2.0
        1.0 3.0 3.0
    ]
    cs_two = [
        0.0 2.0
        2.0 0.0
        1.0 4.0
    ]
    @test_throws AssertionError displacements(cs_one, cs_two)

    disps = displacements(struc_1SSU[5], struc_1SSU[10])
    @test isa(disps, Vector{Float64})
    @test length(disps) == 756
    @test isapprox(disps[20], sqrt(1.984766))
    disps = displacements(struc_1SSU[5], struc_1SSU[10], calphaselector)
    @test length(disps) == 51
    @test isapprox(disps[20], sqrt(0.032822))


    # Test sqdistance and distance
    atom_a = Atom(100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ")
    atom_b = Atom(110, "CA", ' ', [0.0, -1.0, 3.0], 1.0, 10.0, " C", "  ")
    @test sqdistance(atom_a, atom_b) == 10.0
    @test isapprox(distance(atom_a, atom_b), sqrt(10))

    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B']), sqrt(6.852947))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'][50]), sqrt(530.645746))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'][50]["CA"]), sqrt(574.699125))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'], backboneselector), sqrt(17.350083))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'], stdatomselector), sqrt(11.252973))
    @test isapprox(distance(struc_1AKE['A'][50]["CA"], struc_1AKE['B'][50]["CA"]), sqrt(2607.154834))


    # Test bondangle
    atom_a = Atom(100, "CA", ' ', [1.0, 0.0, 1.0], 1.0, 10.0, " C", "  ")
    atom_b = Atom(100, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 10.0, " C", "  ")
    atom_c = Atom(100, "CA", ' ', [3.0, 2.0, 1.0], 1.0, 10.0, " C", "  ")
    @test isapprox(bondangle(atom_a, atom_b, atom_c), 0.713724378944765)
    vec_a = [2.0, 0.0, 0.0]
    vec_b = [2.0, 1.0, 1.0]
    @test isapprox(bondangle(vec_a, vec_b), 0.615479708670387)


    # Test dihedralangle
    atom_a = Atom(100, "CA", ' ', [-1.0, -1.0, 0.0], 1.0, 10.0, " C", "  ")
    atom_b = Atom(100, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 10.0, " C", "  ")
    atom_c = Atom(100, "CA", ' ', [1.0, 0.0, 0.0], 1.0, 10.0, " C", "  ")
    atom_d = Atom(100, "CA", ' ', [2.0, 1.0, -1.0], 1.0, 10.0, " C", "  ")
    @test isapprox(dihedralangle(atom_a, atom_b, atom_c, atom_d), 2.356194490192345)
    vec_a = [1.0, 1.0, 0.0]
    vec_b = [1.0, 0.0, 0.0]
    vec_c = [1.0, -1.0, 1.0]
    @test isapprox(dihedralangle(vec_a, vec_b, vec_c), -0.785398163397448)


    # Test contactmap



end

end # TestStructure
