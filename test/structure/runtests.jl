module TestStructure

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Structure
using Bio.Structure: parsestrict, parselenient, parsevalue, spacestring

const test_files = "test/structure/test_files"


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
    @test resid(disordered_atom) == "10"
    @test resid(disordered_atom, full=true) == "10:A"

    """defaultaltlocid!"""

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
    @test resid(res) == "10"
    @test resid(res, full=true) == "10:A"

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
    @test resid(disordered_res) == "10"
    @test resid(disordered_res, full=true) == "10:A"

    """defaultaltlocid!"""

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

    # Test Chain getters/setters
    @test chainid(chain) == 'A'
    @test resids(chain) == ["10", "H_11"]
    @test sortresids(residues(chain)) == ["10", "H_11"]
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
        1 => Model(1, Dict(
            'A' => Chain('A', ["10"], Dict(
                "10" => Residue("ALA", 'A', 10, ' ', false, ["CA"], Dict(
                    "CA" => Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", "")
                ))
            ))
        ))
    ))

    # Test ProteinStructure getters/setters
    @test structurename(struc) == "test"
    @test modelnumbers(struc) == [1, 5]
    @test isa(models(struc), Dict{Int, Model})
    @test length(models(struc)) == 2
    @test resids(models(struc)[1]['A']) == ["10"]

    """structurename!(struc, "new name")
    @test structurename(struc) == "new name"
    @test modelnumbers(struc) == [1, 5]"""

    # Test ProteinStructure indices
    @test isa(struc[5], Model)
    @test isa(struc[1], Model)
    @test isa(struc['A'], Chain)
    @test_throws KeyError struc[2]
    @test_throws KeyError struc['B']
    @test_throws MethodError struc["A"]
    @test ishetres(struc[5]['B']["H_20"])
    @test !ishetres(struc[1]['A']["10"])
    @test !ishetres(struc['A']["10"])

    # Test ProteinStructure iteration
    struc_list = collect(struc)
    @test isa(struc_list, Vector{Model})
    @test length(struc_list) == 2
    @test modelnumber(struc_list[2]) == 5


    # Test selector functions


    # Test collectatoms and collectresidues


    # Test count functions


    # Test isdisordered


    # Test organise functions



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
    @test_throws PDBParseException parsestrict(line_a, (31,38), Float64, "could not read x coordinate", 10)
    @test_throws PDBParseException parsestrict(line_b, (31,38), Float64, "could not read x coordinate", 10)
    @test_throws PDBParseException parsestrict(line, (7,11), Bool, "could not read atom serial number", 10)

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
    @test_throws PDBParseException parseatomrecord(line_c)
    @test_throws PDBParseException parseatomrecord(line_d)
    @test_throws AssertionError parseatomrecord(line_e)



    struc_1AKE = read("$test_files/1AKE.pdb", PDB)

end


@testset "Writing" begin

end


@testset "Spatial" begin
    # Tests for coordarray



    # Tests for rmsd
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


    # Tests for disps
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

    @test isapprox(disps(coords_one, coords_two), [1.0, sqrt(5)])
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
    @test_throws AssertionError disps(coords_one, coords_two)


    # Tests for dist


end


end # TestStructure
