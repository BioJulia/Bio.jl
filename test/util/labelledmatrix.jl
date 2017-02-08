module TestLSM

using Base.Test

using Bio.Util

@testset "LabelledSquareMatrices" begin

    integral_matrix = """\
    \tA\tB\tC
    A\t1\t2\t4
    B\t2\t1\t2
    C\t4\t2\t1
    """
    integral_matrix_expt = Float64[1 2 4; 2 1 2; 4 2 1;]
    commented_matrix = """\
    # This is a comment that gets skipped. The next blank line should also get
    # skipped.

    \tA\tB\tC
    A\t1\t2\t4
    B\t2\t1\t2
    C\t4\t2\t1
    """
    spacy_matrix = """\
    # This matrix has heaps of extra spaces in the road, to try and trip up the
    # parser.

       A  B  C 
    A 1 2          4
    B  2  1  2   
    C  4  2  1  
    """

    float_matrix = """\
    \tA\tB\tC
    A\t1.0\t2.5\t4.5
    B\t2.5\t1.0\t2.5
    C\t4.5\t2.5\t1.0
    """
    float_matrix_expt = Float64[1. 2.5 4.5; 2.5 1. 2.5; 4.5 2.5 1.;]

    good_matricies = [(float_matrix, float_matrix_expt),
                      (integral_matrix, integral_matrix_expt),
                      (commented_matrix, integral_matrix_expt),
                      (spacy_matrix, integral_matrix_expt),
                     ]

    bad_labels = """\
    \tA\tB\tC
    X\t1\t2\t4
    Y\t2\t1\t2
    Z\t4\t2\t1
    """

    not_square = """\
    \tA\tB\tC
    A\t1\t2\t4
    """

    missing_vals = """\
    \tA\tB\tC
    A\t1\t2\t4
    B\t2\t1\t2
    C\t4
    """
    bad_matricies = [bad_labels,
                     not_square,
                     missing_vals,
                     ]

    file_matricies = [("ints.tab", integral_matrix_expt),
                      ("commented.tab", integral_matrix_expt),
                     ]

    # Maps filename to labels
    NUC44 = ["A", "T", "G", "C", "S", "W", "R", "Y", "K", "M", "B", "V", "H", "D", "N"]
    BLOSUM = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
              "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*"]
    subst_matricies = [("BLOSUM62", BLOSUM),
                       ("NUC.4.4", NUC44),
                     ]
    all_subst_matrices = ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80",
                          "BLOSUM90", "NUC.4.4", "PAM250", "PAM30", "PAM70"]

    @testset "types" for T in [Int64, Float64]
        m, l = readlsm(T, IOBuffer(integral_matrix))
        expt_m = Array{T}(integral_matrix_expt)
        @test m == expt_m
        @test eltype(m) == T
        @test l == ["A", "B", "C"]
    end

    @testset "good_matricies" for (mat, expected) in good_matricies
        m, l = readlsm(IOBuffer(mat))
        @test m == expected
        @test l == ["A", "B", "C"]
    end

    @testset "bad_matricies" for mat in bad_matricies
        @test_throws Exception readlsm(IOBuffer(mat))
    end

    @testset "files" for (fname, expected) in file_matricies
        fname = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "LSM", fname)
        m, l = readlsm(fname)
        @test m == expected
        @test l == ["A", "B", "C"]
    end

    @testset "subst_matrices" for (fname, labels) in subst_matricies
        fname = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "LSM", fname)
        N = length(labels)
        m, l = readlsm(Int64, fname)
        @test size(m) == (N, N)
        @test l == labels
    end

    @testset "all_subst_matrices" for fname in all_subst_matrices
        # Check that we can parse all matrices under the alignm module
        fname = joinpath(dirname(@__FILE__), "..", "..", "src", "align", "data", "submat", fname)
        m, l = readlsm(Int64, fname)
        @test eltype(m) === Int64
    end
end

end
