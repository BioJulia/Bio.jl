module TestStructure

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Structure


@testset "Model" begin
    atom = Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", "")

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
    @test tempfac(atom) == 0.0
    @test element(atom) == "C"
    @test charge(atom) == ""
    @test !ishetero(atom)
    @test resid(atom) == "10"
    @test resid(atom, full=true) == "10:A"

    atom_list = collect(atom)
    @test length(atom_list) == 1
    @test isa(atom_list, Array{Atom,1})
    @test serial(atom_list[1]) == 100


    disordered_atom = DisorderedAtom(Dict(
        'A' => Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.6, 10.0, "C", ""),
        'B' => Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [11.0, 12.0, 13.0], 0.4, 20.0, "C", ""),
    ), 'A')

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

    disordered_atom_list = collect(disordered_atom)
    @test length(disordered_atom_list) == 2
    @test isa(disordered_atom_list, Array{Atom,1})
    @test serial(disordered_atom_list[2]) == 101


    res = Residue("ALA", 'A', 10, ' ', false, ["CA", "CB", "CG"], Dict(
        "CA" => DisorderedAtom(Dict(
            'A' => Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 0.0, "C", ""),
            'B' => Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 0.0, "C", ""),
        ), 'A'),
        "CB" => Atom(false, 102, "CB", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
        "CG" => Atom(false, 103, "CG", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
    ))

    @test resname(res) == "ALA"
    @test chainid(res) == 'A'
    @test resnumber(res) == 10
    @test inscode(res) == ' '
    @test !ishetres(res)
    @test atomnames(res) == ["CA", "CB", "CG"]
    @test resid(res) == "10"
    @test resid(res, full=true) == "10:A"
    @test !ishetero(res)
end


@testset "Parsing" begin

end


@testset "Conversion" begin

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
