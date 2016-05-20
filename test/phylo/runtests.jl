module TestPhylo


if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Phylo
using LightGraphs

@testset "Phylogenies" begin

    @testset "Creation" begin
        # Create an unresolved star phylogeny.
        tips = [:A, :B, :C]
        tree = Phylogeny(tips)

        @testset "Roots" begin
            @test isrooted(tree) == false
            @test isrerootable(tree)
            @test root(tree) == 4
        end

        @testset "Leaves" begin
            @test leaves(tree) == 1:3
            @test isleaf(tree, 1) == true
            @test isleaf(tree, 2) == true
            @test isleaf(tree, 3) == true
        end

        @testset "Clades" begin
            @test isclade(tree, 4) == true
            @test clades(tree) == 5:5
        end

        @testset "Parent - Child" begin
            for i in 1:3
                hasp = i < 4 ? true : false
                nop = i < 4 ? false : true
                nchild = i == 4 ? 3 : 0
                hasc = i == 4 ? true : false
                noc = i == 4 ? false : true
                @test hasparent(tree, i) == hasp
                @test !hasparent(tree, i) == nop
                @test nchildren(tree, i) == nchild
                @test haschildren(tree, i) == hasc
                @test !haschildren(tree, i) == noc
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
        end
        g = DiGraph(5)
        for i in 1:3
            add_edge!(g, Edge(4, i))
        end
        @test tree.graph == g
        @test tree.ntaxa == 3
        @test nv(tree.graph) == 5
    end

end


end # TestPhylo
