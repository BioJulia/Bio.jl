module TestPhylo

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Phylo
using LightGraphs
using Bio.Phylo.Dating
using Bio.Seq

@testset "Phylogenies" begin

    @testset "Creation" begin
        # Create an unresolved star phylogeny.
        tips = [:A, :B, :C]
        stringtips = ["A", "B", "C"]
        tree = Phylogeny(tips)
        tree2 = Phylogeny(stringtips)

        names = [:A, :B, :C, :Root]

        @testset "Roots" begin
            @test isrooted(tree) == false
            @test isrerootable(tree)
            @test root(tree) == 4
            for i in 1:5
                r = i == 4
                @test isroot(tree, i) == r
            end
            for i in names
                r = i == :Root
                @test isroot(tree, i) == r
            end
        end

        @testset "Leaves" begin
            @test leaves(tree) == 1:3
            for i in 1:5
                l = i < 4
                @test isleaf(tree, i) == l
            end
            for i in names
                l = i != :Root
                @test isleaf(tree, i) == l
            end
        end

        @testset "Clades" begin
            @test clades(tree) == 4:5
            for i in 1:5
                c = i > 3
                p = i == 4
                @test isclade(tree, i) == c
                @test ispreterminal(tree, i) == p
                @test issemipreterminal(tree, i) == false
            end
            for i in names
                c = i == :Root
                @test isclade(tree, i) == c
            end
        end

        @testset "Parent - Child" begin
            for i in 1:5
                hasp = i < 4
                nop = i >= 4
                hasp2 = i < 4
                nchild = i == 4 ? 3 : 0
                hasc = i == 4
                hasc2 = i < 4
                noc = i != 4
                @test hasparent(tree, i) == hasp
                @test !hasparent(tree, i) == nop
                @test hasparent(tree, i, 4) == hasp2
                @test hasparent(tree, i, 5) == false
                @test nchildren(tree, i) == nchild
                @test haschildren(tree, i) == hasc
                @test !haschildren(tree, i) == noc
                @test haschild(tree, 4, i) == hasc2
                @test isredundant(tree, i) == false
                if i < 4
                    @test Phylo.parent(tree, i) == 4
                    @test children(tree, i) == Int[]
                else
                    @test_throws ErrorException Phylo.parent(tree, i)
                end
            end
            for i in names
                p = i != :Root
                nchild = i != :Root ? 0 : 3
                @test hasparent(tree, i) == p
                @test nchildren(tree, i) == nchild
            end
            @test hasparent(tree, collect(1:5)) == [true, true, true, false, false]
            @test hasparent(tree) == [true, true, true, false, false]
            @test !hasparent(tree, collect(1:5)) == [false, false, false, true, true]
            @test !hasparent(tree) == [false, false, false, true, true]
            @test nchildren(tree, collect(1:5)) == [0, 0, 0, 3, 0]
            @test nchildren(tree) == [0, 0, 0, 3, 0]
            @test haschildren(tree, collect(1:5)) == [false, false, false, true, false]
            @test haschildren(tree) == [false, false, false, true, false]
            @test !haschildren(tree, collect(1:5)) == [true, true, true, false, true]
            @test !haschildren(tree) == [true, true, true, false, true]
            @test isredundant(tree, collect(1:5)) == [false, false, false, false, false]
            @test isredundant(tree) == [false, false, false, false, false]
        end
        g = DiGraph(5)
        for i in 1:3
            add_edge!(g, Edge(4, i))
        end
        @test tree.graph == g
        @test tree.ntaxa == 3
        @test nv(tree.graph) == 5

        @testset "Copying" begin
            @test copy(tree) == tree
            @test copy(tree) !== tree
        end

        @testset "Branches" begin
            @test branchlength(branchlength!(tree, Edge(4, 1), 0.237), Edge(4, 1)) == 0.237
            @test parent_branch(tree, 1) == Edge(4, 1)
            @test child_branches(tree, 4) == [Edge(4, 1), Edge(4, 2), Edge(4, 3)]
            @test child_branches(rem_branch!(tree, Edge(4, 1)), 4) == [Edge(4, 2), Edge(4, 3)]
            @test child_branches(add_branch!(tree, Edge(4, 1), 0.5), 4) == [Edge(4, 1), Edge(4, 2), Edge(4, 3)]
        end

        @testset "Manipulation" begin
            @test Phylo.unconnected_clades(tree) == [5]
            @test Phylo.subtree_roots(tree) == [4]
            @test child_branches(Phylo.delete!(tree, 1), 4) == [Edge(4, 2), Edge(4, 3)]
            @test child_branches(Phylo.disconnect_root!(tree), 4) == []
        end

    end

end

@testset "Misc Functions" begin
    @test n_possible_rooted(3) == 3
    @test n_possible_rooted(4) == 15
    @test n_possible_rooted(5) == 105
    @test n_possible_unrooted(3) == 1
    @test n_possible_unrooted(4) == 3
    @test n_possible_unrooted(5) == 15
end


@testset "Dating algorithms" begin
    for i in 1:10
        size = rand(10:100)
        mutations = rand(1:(size - 1))
        rate = rand(10e-9:10e-10:10e-6)
        expected = coaltime(size, mutations, rate, SimpleEstimate)
        estimated = coaltime(size, mutations, rate, SpeedDating)
        @test expected in estimated
    end
end

@testset "Distance computation" begin

    # First boolean is test of ANY ambiguous nucleotides.
    # Second boolean is test of same/different ambiguity.
    mutationTests = [
        (DNA_A, DNA_A, false, false),
        (DNA_A, DNA_C, false, false),
        (DNA_A, DNA_G, false, false),
        (DNA_A, DNA_T, false, false),
        (DNA_A, DNA_M, true, true),
        (DNA_A, DNA_R, true, true),
        (DNA_A, DNA_W, true, true),
        (DNA_A, DNA_S, true, false),
        (DNA_A, DNA_Y, true, false),
        (DNA_A, DNA_K, true, false),
        (DNA_A, DNA_V, true, true),
        (DNA_A, DNA_H, true, true),
        (DNA_A, DNA_D, true, true),
        (DNA_A, DNA_B, true, false),
        (DNA_A, DNA_N, true, true),
        (DNA_A, DNA_Gap, true, true),

        (DNA_C, DNA_C, false, false),
        (DNA_C, DNA_G, false, false),
        (DNA_C, DNA_T, false, false),
        (DNA_C, DNA_M, true, true),
        (DNA_C, DNA_R, true, false),
        (DNA_C, DNA_W, true, false),
        (DNA_C, DNA_S, true, true),
        (DNA_C, DNA_Y, true, true),
        (DNA_C, DNA_K, true, false),
        (DNA_C, DNA_V, true, true),
        (DNA_C, DNA_H, true, true),
        (DNA_C, DNA_D, true, false),
        (DNA_C, DNA_B, true, true),
        (DNA_C, DNA_N, true, true),
        (DNA_C, DNA_Gap, true, true),

        (DNA_G, DNA_G, false, false),
        (DNA_G, DNA_T, false, false),
        (DNA_G, DNA_M, true, false),
        (DNA_G, DNA_R, true, true),
        (DNA_G, DNA_W, true, false),
        (DNA_G, DNA_S, true, true),
        (DNA_G, DNA_Y, true, false),
        (DNA_G, DNA_K, true, true),
        (DNA_G, DNA_V, true, true),
        (DNA_G, DNA_H, true, false),
        (DNA_G, DNA_D, true, true),
        (DNA_G, DNA_B, true, true),
        (DNA_G, DNA_N, true, true),
        (DNA_G, DNA_Gap, true, true),

        (DNA_T, DNA_T, false, false),
        (DNA_T, DNA_M, true, false),
        (DNA_T, DNA_R, true, false),
        (DNA_T, DNA_W, true, true),
        (DNA_T, DNA_S, true, false),
        (DNA_T, DNA_Y, true, true),
        (DNA_T, DNA_K, true, true),
        (DNA_T, DNA_V, true, false),
        (DNA_T, DNA_H, true, true),
        (DNA_T, DNA_D, true, true),
        (DNA_T, DNA_B, true, true),
        (DNA_T, DNA_N, true, true),
        (DNA_T, DNA_Gap, true, true),

        (DNA_M, DNA_M, true, true),
        (DNA_M, DNA_R, true, true),
        (DNA_M, DNA_W, true, true),
        (DNA_M, DNA_S, true, true),
        (DNA_M, DNA_Y, true, true),
        (DNA_M, DNA_K, true, false),
        (DNA_M, DNA_V, true, true),
        (DNA_M, DNA_H, true, true),
        (DNA_M, DNA_D, true, true),
        (DNA_M, DNA_B, true, true),
        (DNA_M, DNA_N, true, true),
        (DNA_M, DNA_Gap, true, true),

        (DNA_R, DNA_R, true, true),
        (DNA_R, DNA_W, true, true),
        (DNA_R, DNA_S, true, true),
        (DNA_R, DNA_Y, true, false),
        (DNA_R, DNA_K, true, true),
        (DNA_R, DNA_V, true, true),
        (DNA_R, DNA_H, true, true),
        (DNA_R, DNA_D, true, true),
        (DNA_R, DNA_B, true, true),
        (DNA_R, DNA_N, true, true),
        (DNA_R, DNA_Gap, true, true),

        (DNA_W, DNA_W, true, true),
        (DNA_W, DNA_S, true, false),
        (DNA_W, DNA_Y, true, true),
        (DNA_W, DNA_K, true, true),
        (DNA_W, DNA_V, true, true),
        (DNA_W, DNA_H, true, true),
        (DNA_W, DNA_D, true, true),
        (DNA_W, DNA_B, true, true),
        (DNA_W, DNA_N, true, true),
        (DNA_W, DNA_Gap, true, true),

        (DNA_S, DNA_S, true, true),
        (DNA_S, DNA_Y, true, true),
        (DNA_S, DNA_K, true, true),
        (DNA_S, DNA_V, true, true),
        (DNA_S, DNA_H, true, true),
        (DNA_S, DNA_D, true, true),
        (DNA_S, DNA_B, true, true),
        (DNA_S, DNA_N, true, true),
        (DNA_S, DNA_Gap, true, true),

        (DNA_Y, DNA_Y, true, true),
        (DNA_Y, DNA_K, true, true),
        (DNA_Y, DNA_V, true, true),
        (DNA_Y, DNA_H, true, true),
        (DNA_Y, DNA_D, true, true),
        (DNA_Y, DNA_B, true, true),
        (DNA_Y, DNA_N, true, true),
        (DNA_Y, DNA_Gap, true, true),

        (DNA_K, DNA_K, true, true),
        (DNA_K, DNA_V, true, true),
        (DNA_K, DNA_H, true, true),
        (DNA_K, DNA_D, true, true),
        (DNA_K, DNA_B, true, true),
        (DNA_K, DNA_N, true, true),
        (DNA_K, DNA_Gap, true, true),

        (DNA_V, DNA_V, true, true),
        (DNA_V, DNA_H, true, true),
        (DNA_V, DNA_D, true, true),
        (DNA_V, DNA_B, true, true),
        (DNA_V, DNA_N, true, true),
        (DNA_V, DNA_Gap, true, true),

        (DNA_H, DNA_H, true, true),
        (DNA_H, DNA_D, true, true),
        (DNA_H, DNA_B, true, true),
        (DNA_H, DNA_N, true, true),
        (DNA_H, DNA_Gap, true, true),

        (DNA_D, DNA_D, true, true),
        (DNA_D, DNA_B, true, true),
        (DNA_D, DNA_N, true, true),
        (DNA_D, DNA_Gap, true, true),

        (DNA_B, DNA_B, true, true),
        (DNA_B, DNA_N, true, true),
        (DNA_B, DNA_Gap, true, true),

        (DNA_N, DNA_N, true, true),
        (DNA_N, DNA_Gap, true, true)
    ]

    for test in mutationTests
        @test Phylo.is_ambiguous_strict(test[1], test[2]) == test[3]
        @test Phylo.is_ambiguous_leinient(test[1], test[2], DifferentBases) == test[4]
    end
end

end # Module TestPhylo
