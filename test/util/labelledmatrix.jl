module TestLSM

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

integral_matrix = """\
\tA\tB\tC
A\t1\t2\t4
B\t2\t1\t2
C\t4\t2\t1
"""
integral_matrix_expt = Float64[1 2 4; 2 1 2; 4 2 1;]
commented_matrix = """\
# This is a comment that gets skipped
\tA\tB\tC
A\t1\t2\t4
B\t2\t1\t2
C\t4\t2\t1
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

using Bio.LabelledSquareMatrices

@testset "LabelledSquareMatrices" begin
    @testset "types" for T in [Int64, Float64]
        m, l = readlsm(T, IOBuffer(integral_matrix))
        expt_m = Array{T}(integral_matrix_expt)
        @test m == expt_m
        @test l == UTF8String["A", "B", "C"]
    end

    @testset "good_matricies" for (mat, expected) in good_matricies
        m, l = readlsm(IOBuffer(mat))
        @test m == expected
        @test l == UTF8String["A", "B", "C"]
    end

    @testset "bad_matricies" for mat in bad_matricies
        @test_throws Exception readlsm(IOBuffer(mat))
    end
end

end
