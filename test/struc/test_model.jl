@testset "Model tests" begin
    atom = Atom(false, 100, "CA", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", "")

    disordered_atom = DisorderedAtom(Dict(
        'A' => Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 0.0, "C", ""),
        'B' => Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 0.0, "C", ""),
    ), 'A')

    residue = Residue("ALA", 'A', 10, ' ', false, Dict(
        "CA" => DisorderedAtom(Dict(
            'A' => Atom(false, 100, "CA", 'A', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 0.0, "C", ""),
            'B' => Atom(false, 101, "CA", 'B', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 0.5, 0.0, "C", ""),
        ), 'A'),
        "CB" => Atom(false, 102, "CB", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
        "CG" => Atom(false, 103, "CG", ' ', "ALA", 'A', 10, ' ', [1.0, 2.0, 3.0], 1.0, 0.0, "C", ""),
    ))
end
