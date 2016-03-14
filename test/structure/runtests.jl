module TestStructure

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Structure,
    TestFunctions.get_bio_fmt_specimens
using Bio.Structure: atomid,
    chainidisless,
    parsestrict,
    parselenient,
    parsevalue,
    spacestring


get_bio_fmt_specimens()

# Access a PDB file in BioFmtSpecimens
pdbfilepath(filename::AbstractString) = Pkg.dir("Bio", "test", "BioFmtSpecimens", "PDB", filename)


@testset "Model" begin
    # Test Atom constructor
    atom = Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")

    # Test Atom getters/setters
    @test !ishetatom(atom)
    @test serial(atom) == 100
    @test atomname(atom) == "CA"
    @test altlocid(atom) == ' '
    @test resname(atom) == "ALA"
    @test chainid(atom) == 'A'
    @test resnumber(atom) == 10
    @test inscode(atom) == ' '
    @test x(atom) == 1.0
    @test y(atom) == 2.0
    @test z(atom) == 3.0
    @test coords(atom) == [1.0, 2.0, 3.0]
    @test occupancy(atom) == 1.0
    @test tempfac(atom) == 10.0
    @test element(atom) == "C"
    @test charge(atom) == ""

    @test !ishetero(atom)
    @test !isdisorderedatom(atom)
    @test atomid(atom) == ("10:A", "ALA", "CA")
    @test resid(atom) == "10"
    @test resid(atom, full=true) == "10:A"

    x!(atom, 10.0)
    @test coords(atom) == [10.0, 2.0, 3.0]
    y!(atom, 20.0)
    @test coords(atom) == [10.0, 20.0, 3.0]
    z!(atom, 30.0)
    @test coords(atom) == [10.0, 20.0, 30.0]
    coords!(atom, [40.0, 50.0, 60.0])
    @test coords(atom) == [40.0, 50.0, 60.0]

    # Test Atom iteration
    atom_list = collect(atom)
    @test isa(atom_list, Vector{Atom})
    @test length(atom_list) == 1
    @test serial(atom_list[1]) == 100

    # Test DisorderedAtom constructor
    disordered_atom = DisorderedAtom(Dict(
        'A' => Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.6, 10.0, "C", ""),
        'B' => Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [11.0, 12.0, 13.0], 0.4, 20.0, "C", "")
    ), 'A')

    # Test DisorderedAtom getters/setters
    @test defaultaltlocid(disordered_atom) == 'A'
    @test isa(defaultatom(disordered_atom), Atom)
    @test serial(defaultatom(disordered_atom)) == 100
    @test altlocids(disordered_atom) == ['A', 'B']

    @test !ishetatom(disordered_atom)
    @test serial(disordered_atom) == 100
    @test atomname(disordered_atom) == "CA"
    @test altlocid(disordered_atom) == 'A'
    @test resname(disordered_atom) == "ALA"
    @test chainid(disordered_atom) == 'A'
    @test resnumber(disordered_atom) == 10
    @test inscode(disordered_atom) == ' '
    @test x(disordered_atom) == 1.0
    @test y(disordered_atom) == 2.0
    @test z(disordered_atom) == 3.0
    @test coords(disordered_atom) == [1.0, 2.0, 3.0]
    @test occupancy(disordered_atom) == 0.6
    @test tempfac(disordered_atom) == 10.0
    @test element(disordered_atom) == "C"
    @test charge(disordered_atom) == ""

    @test !ishetero(disordered_atom)
    @test isdisorderedatom(disordered_atom)
    @test atomid(disordered_atom) == ("10:A", "ALA", "CA")
    @test resid(disordered_atom) == "10"
    @test resid(disordered_atom, full=true) == "10:A"

    disordered_atom_mod = DisorderedAtom(disordered_atom, 'B')
    @test defaultaltlocid(disordered_atom_mod) == 'B'
    @test serial(disordered_atom_mod) == 101
    @test_throws AssertionError DisorderedAtom(disordered_atom, 'C')

    # Only the coordinates on the default atom are changed
    x!(disordered_atom, 10.0)
    @test coords(disordered_atom) == [10.0, 2.0, 3.0]
    @test coords(disordered_atom['B']) == [11.0, 12.0, 13.0]
    y!(disordered_atom, 20.0)
    @test coords(disordered_atom) == [10.0, 20.0, 3.0]
    @test coords(disordered_atom['B']) == [11.0, 12.0, 13.0]
    z!(disordered_atom, 30.0)
    @test coords(disordered_atom) == [10.0, 20.0, 30.0]
    @test coords(disordered_atom['B']) == [11.0, 12.0, 13.0]
    coords!(disordered_atom, [40.0, 50.0, 60.0])
    @test coords(disordered_atom) == [40.0, 50.0, 60.0]
    @test coords(disordered_atom['B']) == [11.0, 12.0, 13.0]
    x!(disordered_atom['A'], 100.0)
    @test coords(disordered_atom) == [100.0, 50.0, 60.0]
    @test coords(disordered_atom['B']) == [11.0, 12.0, 13.0]
    x!(disordered_atom['B'], 110.0)
    @test coords(disordered_atom) == [100.0, 50.0, 60.0]
    @test coords(disordered_atom['B']) == [110.0, 12.0, 13.0]

    # Test DisorderedAtom indices
    @test isa(disordered_atom['A'], Atom)
    @test serial(disordered_atom['A']) == 100
    @test serial(disordered_atom['B']) == 101
    @test_throws KeyError disordered_atom['C']
    @test_throws MethodError disordered_atom["A"]

    # Test DisorderedAtom iteration
    disordered_atom_list = collect(disordered_atom)
    @test isa(disordered_atom_list, Vector{Atom})
    @test length(disordered_atom_list) == 2
    @test serial(disordered_atom_list[2]) == 101


    # Test Residue constructor
    res = Residue("ALA", 'A', 10, ' ', false, ["CA", "CB", "CG"], Dict(
        "CA" => disordered_atom,
        "CB" => Atom(false, 102, "CB", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
        "CG" => Atom(false, 103, "CG", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", "")
    ))
    res_min = Residue("ALA", 'A', 10, ' ', false)
    res_min = Residue(AbstractAtom[
        Atom(false, 102, "CG", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
        Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
        Atom(false, 101, "CB", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", "")
    ])

    # Test Residue getters/setters
    @test resname(res) == "ALA"
    @test chainid(res) == 'A'
    @test resnumber(res) == 10
    @test inscode(res) == ' '
    @test !ishetres(res)
    @test atomnames(res) == ["CA", "CB", "CG"]
    @test isa(atoms(res), Dict{ASCIIString, AbstractAtom})
    @test length(atoms(res)) == 3
    @test serial(atoms(res)["CA"]) == 100

    @test !ishetero(res)
    @test !isdisorderedres(res)
    @test resid(res) == "10"
    @test resid(res, full=true) == "10:A"

    @test resname(res_min) == "ALA"
    @test chainid(res_min) == 'A'
    @test resnumber(res_min) == 10
    @test inscode(res_min) == ' '
    @test !ishetres(res_min)
    @test atomnames(res_min) == ["CA", "CB", "CG"]
    @test isa(atoms(res_min), Dict{ASCIIString, AbstractAtom})
    @test length(atoms(res_min)) == 3
    @test serial(atoms(res_min)["CA"]) == 100

    @test !ishetero(res_min)
    @test !isdisorderedres(res_min)
    @test resid(res_min) == "10"
    @test resid(res_min, full=true) == "10:A"

    # Test Residue indices
    @test isa(res["CA"], AbstractAtom)
    @test serial(res["CA"]) == 100
    @test serial(res["CB"]) == 102
    @test_throws KeyError res["N"]

    # Test Residue iteration
    res_list = collect(res)
    @test isa(res_list, Vector{AbstractAtom})
    @test length(res_list) == 3
    @test serial(res_list[3]) == 103

    # Test DisorderedResidue constructor
    disordered_res = DisorderedResidue(Dict(
        "ALA" => res,
        "VAL" => Residue("VAL", 'A', 10, ' ', false, ["CA"], Dict(
            "CA" => Atom(false, 110, "CA", ' ', "VAL", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", "")
        ))
    ), "ALA")

    # Test DisorderedResidue getters/setters
    @test isa(disorderedres(disordered_res, "VAL"), AbstractResidue)
    @test resname(disorderedres(disordered_res, "VAL")) == "VAL"
    @test defaultresname(disordered_res) == "ALA"
    @test isa(defaultresidue(disordered_res), Residue)
    @test resname(defaultresidue(disordered_res)) == "ALA"
    @test resnames(disordered_res) == ["ALA", "VAL"]

    @test resname(disordered_res) == "ALA"
    @test chainid(disordered_res) == 'A'
    @test resnumber(disordered_res) == 10
    @test inscode(disordered_res) == ' '
    @test !ishetres(disordered_res)
    @test atomnames(disordered_res) == ["CA", "CB", "CG"]
    @test isa(atoms(res), Dict{ASCIIString, AbstractAtom})
    @test length(atoms(res)) == 3
    @test serial(atoms(res)["CA"]) == 100

    @test !ishetero(disordered_res)
    @test isdisorderedres(disordered_res)
    @test resid(disordered_res) == "10"
    @test resid(disordered_res, full=true) == "10:A"

    disordered_res_mod = DisorderedResidue(disordered_res, "VAL")
    @test defaultresname(disordered_res_mod) == "VAL"
    @test atomnames(disordered_res_mod) == ["CA"]
    @test_throws AssertionError DisorderedResidue(disordered_res, "SER")

    # Test DisorderedResidue indices
    @test isa(disordered_res["CA"], AbstractAtom)
    @test serial(disordered_res["CA"]) == 100
    @test serial(disordered_res["CB"]) == 102
    @test_throws KeyError res["N"]

    # Test DisorderedResidue iteration
    disordered_res_list = collect(disordered_res)
    @test isa(res_list, Vector{AbstractAtom})
    @test length(disordered_res_list) == 3
    @test serial(res_list[3]) == 103


    # Test Chain constructor
    chain = Chain('A', ["10", "H_11"], Dict(
        "10" => disordered_res,
        "H_11" => Residue("ATP", 'A', 11, ' ', true, ["PB"], Dict(
            "PB" => Atom(true, 120, "PB", ' ', "ATP", 'A', 11, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "P", "")
        ))
    ))
    chain_min = Chain('A')

    # Test Chain getters/setters
    @test chainid(chain) == 'A'
    @test resids(chain) == ["10", "H_11"]
    @test isa(residues(chain), Dict{ASCIIString, AbstractResidue})
    @test length(residues(chain)) == 2
    @test serial(residues(chain)["10"]["CA"]) == 100

    # Test Chain indices
    @test isa(chain["10"], AbstractResidue)
    @test isa(chain[10], AbstractResidue)
    @test_throws KeyError chain["H_10"]
    @test serial(chain["10"]["CA"]) == 100
    @test serial(chain[10]["CA"]) == 100
    @test resname(chain["H_11"]) == "ATP"
    @test_throws KeyError chain["11"]
    @test_throws KeyError chain[11]

    # Test Chain iteration
    chain_list = collect(chain)
    @test isa(chain_list, Vector{AbstractResidue})
    @test length(chain_list) == 2
    @test resname(chain_list[2]) == "ATP"


    # Test Model constructor
    model = Model(5, Dict(
        'A' => chain,
        'B' => Chain('B', ["H_20"], Dict(
            "H_20" => Residue("ATP", 'B', 1, ' ', true, ["PB"], Dict(
                "PB" => Atom(true, 1000, "PB", ' ', "ATP", 'B', 1, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "P", "")
            ))
        ))
    ))
    model_min = Model(1)
    model_min = Model()

    # Test Model getters/setters
    @test modelnumber(model) == 5
    @test chainids(model) == ['A', 'B']
    @test isa(chains(model), Dict{Char, Chain})
    @test length(chains(model)) == 2
    @test resname(chains(model)['B']["H_20"]) == "ATP"

    # Test Model indices
    @test isa(model['A'], Chain)
    @test resname(model['A']["H_11"]) == "ATP"
    @test_throws KeyError model['C']
    @test_throws MethodError model["A"]

    # Test Model iteration
    model_list = collect(model)
    @test isa(model_list, Vector{Chain})
    @test length(model_list) == 2
    @test chainid(model_list[2]) == 'B'


    # Test ProteinStructure constructor
    struc = ProteinStructure("test", Dict(
        5 => model,
        3 => Model(3, Dict(
            'A' => Chain('A', ["10"], Dict(
                "10" => Residue("ALA", 'A', 10, ' ', false, ["CA"], Dict(
                    "CA" => Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", "")
                ))
            ))
        ))
    ))
    struc_min = ProteinStructure("test")
    struc_min = ProteinStructure()

    # Test ProteinStructure getters/setters
    @test structurename(struc) == "test"
    @test modelnumbers(struc) == [3, 5]
    @test isa(models(struc), Dict{Int, Model})
    @test length(models(struc)) == 2
    @test resids(models(struc)[3]['A']) == ["10"]

    @test isa(defaultmodel(struc), Model)
    @test modelnumber(defaultmodel(struc)) == 3
    @test chainids(struc) == ['A']

    # Test ProteinStructure indices
    @test isa(struc[5], Model)
    @test isa(struc[3], Model)
    @test isa(struc['A'], Chain)
    @test_throws KeyError struc[2]
    @test_throws KeyError struc['B']
    @test_throws MethodError struc["A"]
    @test ishetres(struc[5]['B']["H_20"])
    @test !ishetres(struc[3]['A']["10"])
    @test !ishetres(struc['A']["10"])

    # Test ProteinStructure iteration
    struc_list = collect(struc)
    @test isa(struc_list, Vector{Model})
    @test length(struc_list) == 2
    @test modelnumber(struc_list[2]) == 5


    # Test selector functions
    atom_a = Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    atom_b = Atom(true, 110, "MG", ' ', "MG", 'A', 11, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    res_a = Residue("ALA", 'A', 10, ' ', false, [], Dict())
    res_b = Residue("HOH", 'A', 100, ' ', true, [], Dict())

    @test stdatomselector(atom_a)
    @test !stdatomselector(atom_b)
    @test !hetatomselector(atom_a)
    @test hetatomselector(atom_b)
    @test atomnameselector(atom_a, Set(["CA", "N", "C"]))
    @test atomnameselector(atom_a, ["CA", "N", "C"])
    @test !atomnameselector(atom_b, Set(["CA", "N", "C"]))
    @test calphaselector(atom_a)
    @test !calphaselector(atom_b)
    @test !calphaselector(Atom(true, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", ""))
    @test backboneselector(atom_a)
    @test !backboneselector(atom_b)
    @test !backboneselector(Atom(true, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", ""))
    @test heavyatomselector(atom_a)
    @test !heavyatomselector(atom_b)
    @test !heavyatomselector(Atom(false, 100, "H1", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "H", ""))
    @test resnameselector(atom_a, Set(["ALA"]))
    @test resnameselector(atom_a, ["ALA"])
    @test !resnameselector(atom_b, Set(["ALA"]))
    @test resnameselector(res_a, Set(["ALA"]))
    @test !resnameselector(res_b, Set(["ALA"]))
    @test !waterselector(res_a)
    @test waterselector(res_b)
    @test stdresselector(res_a)
    @test !stdresselector(res_b)
    @test !hetresselector(res_a)
    @test hetresselector(res_b)
    @test !disorderselector(atom_a)
    @test disorderselector(disordered_atom)
    @test !disorderselector(res_a)
    @test disorderselector(disordered_res)
    @test hydrogenselector(Atom(false, 100, "H", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "H", ""))
    @test !hydrogenselector(Atom(false, 100, "H", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", ""))
    @test hydrogenselector(Atom(false, 100, "H1", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "", ""))
    @test hydrogenselector(Atom(false, 100, "1H", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "", ""))
    @test !hydrogenselector(Atom(false, 100, "NH1", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "", ""))


    # Further tests for structural element ordering
    # Order when looping over a DisorderedAtom is the atom serial
    disordered_atom_ord = DisorderedAtom(Dict(
        'A' => Atom(false, 102, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.3, 10.0, "C", ""),
        'B' => Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.4, 10.0, "C", ""),
        'C' => Atom(false, 100, "CA", 'C', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.3, 10.0, "C", ""),
    ), 'B')
    @test altlocids(disordered_atom_ord) == ['C', 'B', 'A']

    # Order when sorting an atom list is the atom serial
    atom_list_ord = AbstractAtom[
        DisorderedAtom(Dict(
            'A' => Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.4, 10.0, "C", ""),
            'B' => Atom(false, 104, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.6, 10.0, "C", "")
        ), 'B'),
        Atom(false, 102, "CB", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", ""),
        Atom(false, 103, "CG", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    ]
    @test map(atomname, sort(atom_list_ord)) == ["CB", "CG", "CA"]
    @test map(atomname, sort(atom_list_ord; rev=true)) == ["CA", "CG", "CB"]
    sort!(atom_list_ord)
    @test map(atomname, atom_list_ord) == ["CB", "CG", "CA"]

    # Order when sorting a residue list is chain ID, then stdres/hetres, then residue number, then insertion code
    residues_ord = AbstractResidue[
        Residue("ALA", 'A', 201, 'A', false),
        Residue("ALA", 'A', 203, ' ', false),
        Residue("ALA", 'A', 200, ' ', true),
        Residue("ALA", 'A', 201, 'B', false),
        Residue("ALA", 'A', 202, ' ', false),
        Residue("ALA", 'B', 300, ' ', false),
        Residue("ALA", 'A', 201, ' ', true),
        Residue("ALA", 'A', 201, ' ', false),
        Residue("ALA", 'A', 201, 'A', true),
        Residue("ALA", 'B', 100, ' ', false),
        Residue("ALA", 'A', 203, ' ', true),
        Residue("ALA", 'A', 200, ' ', false),
    ]
    # Test resid
    @test map(resid, residues_ord) == ["201A", "203", "H_200", "201B", "202", "300", "H_201", "201", "H_201A", "100", "H_203", "200"]
    @test map(res -> resid(res; full=true), residues_ord) == ["201A:A", "203:A", "H_200:A", "201B:A", "202:A", "300:B", "H_201:A", "201:A", "H_201A:A", "100:B", "H_203:A", "200:A"]
    @test map(res -> resid(res; full=true), sort(residues_ord)) == ["200:A", "201:A", "201A:A", "201B:A", "202:A", "203:A", "H_200:A", "H_201:A", "H_201A:A", "H_203:A", "100:B", "300:B"]
    @test map(res -> resid(res; full=true), sort(residues_ord; rev=true)) == ["300:B", "100:B", "H_203:A", "H_201A:A", "H_201:A", "H_200:A", "203:A", "202:A", "201B:A", "201A:A", "201:A", "200:A"]
    sort!(residues_ord)
    @test map(res -> resid(res; full=true), residues_ord) == ["200:A", "201:A", "201A:A", "201B:A", "202:A", "203:A", "H_200:A", "H_201:A", "H_201A:A", "H_203:A", "100:B", "300:B"]

    # Order of listing residue names in a DisorderedResidue is default then alphabetical
    disordered_res_ord = DisorderedResidue(Dict(
        "THR" => Residue("THR", 'A', 201, ' ', false),
        "ALA" => Residue("ALA", 'A', 201, ' ', false),
        "ILE" => Residue("ILE", 'A', 201, ' ', false),
        "SER" => Residue("SER", 'A', 201, ' ', false),
        "VAL" => Residue("VAL", 'A', 201, ' ', false)
    ), "SER")
    @test defaultresname(disordered_res_ord) == "SER"
    @test resnames(disordered_res_ord) == ["SER", "ALA", "ILE", "THR", "VAL"]

    # Order when sorting chain IDs is character ordering with the empty chain ID at the end
    @test chainidisless('A', 'B')
    @test chainidisless('A', 'a')
    @test chainidisless('1', 'A')
    @test chainidisless('A', ' ')
    @test !chainidisless('B', 'A')
    @test !chainidisless(' ', 'A')
    @test !chainidisless(' ', ' ')
    model_ord = Model(1, Dict(
        'A' => Chain('A'),
        ' ' => Chain(' '),
        '1' => Chain('1'),
        'a' => Chain('a'),
        'X' => Chain('X'),
    ))
    @test chainids(model_ord) == ['1', 'A', 'X', 'a', ' ']
    @test sortchainids(['A', ' ', '1', 'a', 'X']) == ['1', 'A', 'X', 'a', ' ']
end


@testset "Parsing" begin
    # Test parsevalue
    line = "ATOM     40  CB  LEU A   5      22.088  45.547  29.675  1.00 22.23           C  "
    @test parsevalue(line, (7,11), Int) == 40
    @test parsevalue(line, (31,38), Float64) == 22.088
    @test strip(parsevalue(line, (13,16), ASCIIString)) == "CB"
    @test parsevalue(line, (22,22), Char) == 'A'
    @test_throws ErrorException parsevalue(line, (1,4), Int)
    @test_throws ErrorException parsevalue(line, (1,4), Bool)
    @test_throws ErrorException parsevalue(line, (79,100), Int)

    # Test parsestrict
    line =   "ATOM    591  C   GLY A  80      29.876  54.131  35.806  1.00 40.97           C  "
    line_a = "ATOM    591  C   GLY A  80              54.131  35.806  1.00 40.97           C  "
    line_b = "ATOM    591  C   GLY A  80 "
    @test parsestrict(line, (7,11), Int, "could not read atom serial number", 10) == 591
    @test strip(parsestrict(line, (13,16), ASCIIString, "could not read atom name", 20)) == "C"
    @test_throws PDBParseError parsestrict(line_a, (31,38), Float64, "could not read x coordinate", 10)
    @test_throws PDBParseError parsestrict(line_b, (31,38), Float64, "could not read x coordinate", 10)
    @test_throws PDBParseError parsestrict(line, (7,11), Bool, "could not read atom serial number", 10)

    # Test parselenient
    line =   "ATOM     40  CB  LEU A   5      22.088  45.547  29.675  1.00 22.23           C  "
    line_a = "ATOM     40  CB  LEU A   5      22.088  45.547  29.675  1.00 22.23              "
    line_b = "ATOM     40  CB  LEU A   5      22.088  45.547  29.675  1.00 22.23  "
    line_c = "ATOM     40  CB  LEU A   5      22.088  45.547  29.675       22.23           C  "
    @test strip(parselenient(line, (77,78), ASCIIString, "")) == "C"
    @test parselenient(line_a, (77,78), ASCIIString, "") == "  "
    @test parselenient(line_b, (77,78), ASCIIString, "") == ""
    @test parselenient(line_b, (77,78), ASCIIString, "C") == "C"
    @test parselenient(line_c, (55,60), Float64, 1.0) == 1.0
    @test parselenient(line, (77,78), Bool, "N") == "N"

    # Test parseatomrecord
    line_a = "ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  "
    line_b = "HETATM 3474  O  B XX A 334A      8.802  62.000   8.672  1.00 39.15           O1-"
    line_c = "ATOM    669  CA  ILE A  90      xxxxxx  33.110  31.221  1.00 25.76           C  "
    line_d = "ATOM    669  CA  ILE A  90      31.743   "
    line_e = "REMARK   1 REFERENCE 1                                                          "
    atom = parseatomrecord(line_a, 10)
    @test !ishetatom(atom)
    @test serial(atom) == 669
    @test atomname(atom) == "CA"
    @test altlocid(atom) == ' '
    @test resname(atom) == "ILE"
    @test chainid(atom) == 'A'
    @test resnumber(atom) == 90
    @test inscode(atom) == ' '
    @test x(atom) == 31.743
    @test y(atom) == 33.110
    @test z(atom) == 31.221
    @test coords(atom) == [31.743, 33.110, 31.221]
    @test occupancy(atom) == 1.00
    @test tempfac(atom) == 25.76
    @test element(atom) == "C"
    @test charge(atom) == ""
    atom = parseatomrecord(line_b)
    @test ishetatom(atom)
    @test serial(atom) == 3474
    @test atomname(atom) == "O"
    @test altlocid(atom) == 'B'
    @test resname(atom) == "XX"
    @test chainid(atom) == 'A'
    @test resnumber(atom) == 334
    @test inscode(atom) == 'A'
    @test x(atom) == 8.802
    @test y(atom) == 62.0
    @test z(atom) == 8.672
    @test occupancy(atom) == 1.00
    @test tempfac(atom) == 39.15
    @test element(atom) == "O"
    @test charge(atom) == "1-"
    @test_throws PDBParseError parseatomrecord(line_c)
    @test_throws PDBParseError parseatomrecord(line_d)
    @test_throws AssertionError parseatomrecord(line_e)

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

    # Test collectatoms
    atoms = collectatoms(struc)
    @test length(atoms) == 3804
    @test isa(atoms, Vector{AbstractAtom})
    @test isa(atoms[70], Atom)
    @test isa(atoms[1290], DisorderedAtom)
    @test serial(atoms[1660]) == 3323
    atoms = collectatoms(struc, hetatomselector)
    @test length(atoms) == 492
    @test serial(atoms[60]) == 3443
    atoms = collectatoms(struc, disorderselector)
    @test length(atoms) == 12
    @test serial(atoms[10]) == 3338
    atoms = collectatoms(struc, stdatomselector, disorderselector)
    @test length(atoms) == 5
    @test serial(atoms[4]) == 1294
    atoms = collectatoms(struc[1])
    @test length(atoms) == 3804
    @test serial(atoms[1660]) == 3323
    atoms = collectatoms(struc['A'])
    @test length(atoms) == 1954
    @test serial(atoms[240]) == 240
    atoms = collectatoms(struc['A'][50])
    @test length(atoms) == 9
    @test serial(atoms[4]) == 358
    atoms = collectatoms(struc['A'][50]["CA"])
    @test length(atoms) == 1
    @test isa(atoms, Vector{AbstractAtom})
    @test serial(atoms[1]) == 356
    atoms = collectatoms(struc['A'][167]["CZ"])
    @test length(atoms) == 1
    @test isa(atoms, Vector{AbstractAtom})
    @test isa(atoms[1], DisorderedAtom)
    @test serial(atoms[1]) == 1292
    atoms = collectatoms(Chain[struc['B'], struc['A']])
    @test length(atoms) == 3804
    @test serial(atoms[5]) == 1667
    atoms = collectatoms(Residue[struc['A'][51], struc['A'][50]])
    @test length(atoms) == 17
    @test serial(atoms[10]) == 356
    atoms = collectatoms(AbstractResidue[struc['A'][51], struc['A'][50]])
    @test length(atoms) == 17
    @test serial(atoms[10]) == 356
    atoms = collectatoms(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test length(atoms) == 2
    @test isa(atoms, Vector{Atom})
    @test serial(atoms[1]) == 356
    atoms = collectatoms(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(atoms) == 2
    @test isa(atoms, Vector{DisorderedAtom})
    @test serial(atoms[2]) == 1292
    atoms = collectatoms(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(atoms) == 2
    @test isa(atoms, Vector{AbstractAtom})
    @test serial(atoms[2]) == 1292


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


    # Test countresidues
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


    # Test formatomlist
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


    # Test parsing 1SSU (multiple models)
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
        collect(struc['A'][9]);
        collect(disorderedres(struc['A'][10], "SER"));
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


    # Test parser error handling
    error = PDBParseError("message", 10, "line")
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
    @test spaceatomname(Atom(false, 1, "N", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "N", "")) == " N  "
    @test spaceatomname(Atom(false, 1, "N", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "", "")) == " N  "
    @test spaceatomname(Atom(false, 1, "CA", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "C", "")) == " CA "
    @test spaceatomname(Atom(false, 1, "NE1", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "N", "")) == " NE1"
    @test spaceatomname(Atom(false, 1, "2HD1", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "H", "")) == "2HD1"
    @test spaceatomname(Atom(false, 1, "HH11", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "H", "")) == "HH11"
    @test spaceatomname(Atom(false, 1, "1H", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "H", "")) == "1H  "
    @test spaceatomname(Atom(false, 1, "MG", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "MG", "")) == "MG  "
    @test spaceatomname(Atom(false, 1, "MG", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "", "")) == " MG "
    @test_throws AssertionError spaceatomname(Atom(false, 1, "11H", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "H", ""))
    @test_throws AssertionError spaceatomname(Atom(false, 1, "11H11", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "H", ""))
    @test_throws AssertionError spaceatomname(Atom(false, 1, "1MG", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "MG", ""))


    # Test pdbline
    # These tests should be changed long term to require the conventional decimal formatting, e.g. 0.50 not 0.5 for occupancy
    atom = Atom(false, 10, "N", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "N", "")
    @test join(pdbline(atom)) == "ATOM     10  N   ALA A   1         0.0     0.0     0.0   1.0   0.0           N  "
    atom = Atom(true, 101, "C", 'A', "LEU", 'B', 20, ' ', [10.5, 20.12345, -5.1227], 0.50, 50.126, "C", "1+")
    line = join(pdbline(atom))
    @test line == "HETATM  101  C  ALEU B  20        10.5  20.123  -5.123   0.5 50.13           C1+"
    atom = parseatomrecord(line)
    @test ishetatom(atom)
    @test serial(atom) == 101
    @test atomname(atom) == "C"
    @test altlocid(atom) == 'A'
    @test resname(atom) == "LEU"
    @test chainid(atom) == 'B'
    @test resnumber(atom) == 20
    @test inscode(atom) == ' '
    @test coords(atom) == [10.5, 20.123, -5.123]
    @test occupancy(atom) == 0.5
    @test tempfac(atom) == 50.13
    @test element(atom) == "C"
    @test charge(atom) == "1+"
    @test_throws AssertionError pdbline(Atom(false, 1, "11H11", ' ', "ALA", 'A', 1, ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "H", ""))


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

    rm(temp_filename)
end


@testset "Spatial" begin
    # Test coordarray
    atom = Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    coords = coordarray(atom)
    @test size(coords) == (3,1)
    @test coords[1] == 1.0
    @test coords[2] == 2.0
    @test coords[3] == 3.0

    struc_1AKE = read(pdbfilepath("1AKE.pdb"), PDB)
    coords = coordarray(struc_1AKE)
    @test size(coords) == (3,3804)
    @test coords[1,3787] == 20.135
    @test coords[2,3787] == -10.789
    @test coords[3,3787] == -1.732
    coords = coordarray(struc_1AKE['A'], calphaselector)
    @test size(coords) == (3,214)
    @test coords[1,10] == 17.487
    @test coords[2,10] == 42.426
    @test coords[3,10] == 19.756
    @test coordarray(coords) == coords


    # Test rmsd
    coords_one = [
        0.0 0.0;
        1.0 0.0;
        1.0 3.0;
    ]
    coords_two = [
        0.0 2.0;
        2.0 0.0;
        1.0 3.0;
    ]
    @test isapprox(rmsd(coords_one, coords_two), sqrt(5/2))
    coords_one = [
        0.0 0.0 1.0;
        1.0 0.0 2.0;
        1.0 3.0 3.0;
    ]
    coords_two = [
        0.0 2.0;
        2.0 0.0;
        1.0 3.0;
    ]
    @test_throws AssertionError rmsd(coords_one, coords_two)

    struc_1SSU = read(pdbfilepath("1SSU.pdb"), PDB)
    @test isapprox(rmsd(struc_1SSU[1], struc_1SSU[2], calphaselector), 4.1821925809691889)
    @test isapprox(rmsd(struc_1SSU[5], struc_1SSU[6], backboneselector), 5.2878196391279939)
    @test_throws AssertionError rmsd(struc_1SSU[1]['A'][8], struc_1SSU[1]['A'][9])


    # Test displacements
    coords_one = [
        0.0 0.0;
        1.0 0.0;
        1.0 3.0;
    ]
    coords_two = [
        0.0 2.0;
        2.0 0.0;
        1.0 4.0;
    ]
    @test isapprox(displacements(coords_one, coords_two), [1.0, sqrt(5)])
    coords_one = [
        0.0 0.0 1.0;
        1.0 0.0 2.0;
        1.0 3.0 3.0;
    ]
    coords_two = [
        0.0 2.0;
        2.0 0.0;
        1.0 4.0;
    ]
    @test_throws AssertionError displacements(coords_one, coords_two)

    disps = displacements(struc_1SSU[5], struc_1SSU[10])
    @test isa(disps, Vector{Float64})
    @test length(disps) == 756
    @test isapprox(disps[20], sqrt(1.984766))
    disps = displacements(struc_1SSU[5], struc_1SSU[10], calphaselector)
    @test length(disps) == 51
    @test isapprox(disps[20], sqrt(0.032822))


    # Test distance
    atom_a = Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "")
    atom_b = Atom(false, 110, "CA", ' ', "ALA", 'A', 11, ' ', [0.0, -1.0, 3.0], 1.0, 10.0, "C", "")
    @test isapprox(distance(atom_a, atom_b), sqrt(10))

    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B']), sqrt(6.852947))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'][50]), sqrt(530.645746))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'][50]["CA"]), sqrt(574.699125))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'], backboneselector), sqrt(17.350083))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'], stdatomselector), sqrt(11.252973))
    @test isapprox(distance(struc_1AKE['A'][50]["CA"], struc_1AKE['B'][50]["CA"]), sqrt(2607.154834))
end

end # TestStructure
