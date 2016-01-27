@testset "Model tests" begin
    atom = Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", "")

    @test !ishetatom(atom)
    @test getserial(atom) == 100
    @test getatomname(atom) == "CA"
    @test getaltlocid(atom) == ' '
    @test getresname(atom) == "ALA"
    @test getchainid(atom) == 'A'
    @test getresnumber(atom) == 10
    @test getinscode(atom) == ' '
    @test getx(atom) == 1.0
    @test gety(atom) == 2.0
    @test getz(atom) == 3.0
    @test getcoords(atom) == [1.0, 2.0, 3.0]
    @test getoccupancy(atom) == 1.0
    @test gettempfac(atom) == 0.0
    @test getelement(atom) == "C"
    @test getcharge(atom) == ""
    @test !ishetero(atom)
    @test getresid(atom) == "10"
    @test getresid(atom, full=true) == "10:A"

    atom_list = collect(atom)
    @test length(atom_list) == 1
    @test typeof(atom_list) == Array{Atom,1}
    @test getserial(atom_list[1]) == 100


    disordered_atom = DisorderedAtom(Dict(
        'A' => Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.6, 10.0, "C", ""),
        'B' => Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [11.0, 12.0, 13.0], 0.4, 20.0, "C", ""),
    ), 'A')

    @test getdefaultaltlocid(disordered_atom) == 'A'
    @test typeof(getdefaultatom(disordered_atom)) == Atom
    @test getserial(getdefaultatom(disordered_atom)) == 100
    @test getaltlocids(disordered_atom) == ['A', 'B']

    @test !ishetatom(disordered_atom)
    @test getserial(disordered_atom) == 100
    @test getatomname(disordered_atom) == "CA"
    @test getaltlocid(disordered_atom) == 'A'
    @test getresname(disordered_atom) == "ALA"
    @test getchainid(disordered_atom) == 'A'
    @test getresnumber(disordered_atom) == 10
    @test getinscode(disordered_atom) == ' '
    @test getx(disordered_atom) == 1.0
    @test gety(disordered_atom) == 2.0
    @test getz(disordered_atom) == 3.0
    @test getcoords(disordered_atom) == [1.0, 2.0, 3.0]
    @test getoccupancy(disordered_atom) == 0.6
    @test gettempfac(disordered_atom) == 10.0
    @test getelement(disordered_atom) == "C"
    @test getcharge(disordered_atom) == ""
    @test !ishetero(disordered_atom)
    @test getresid(disordered_atom) == "10"
    @test getresid(disordered_atom, full=true) == "10:A"

    disordered_atom_list = collect(disordered_atom)
    @test length(disordered_atom_list) == 2
    @test typeof(disordered_atom_list) == Array{Atom,1}
    @test getserial(atom_list[2]) == 101


    res = Residue("ALA", 'A', 10, ' ', false, Dict(
        "CA" => DisorderedAtom(Dict(
            'A' => Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 0.0, "C", ""),
            'B' => Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 0.0, "C", ""),
        ), 'A'),
        "CB" => Atom(false, 102, "CB", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
        "CG" => Atom(false, 103, "CG", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
    ))

    @test getresname(res) == "ALA"
    @test getchainid(res) == 'A'
    @test getresnumber(res) == 10
    @test getinscode(res) == ' '
    @test !ishetres(res)
    @test getatomnames(res) == ["CA", "CB", "CG"]
    @test getresid(res) == "10"
    @test getresid(res, full=true) == "10:A"
    @test !ishetero(res)
end
