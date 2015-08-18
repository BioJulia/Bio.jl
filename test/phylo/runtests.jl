module TestPhylo

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Phylo

@testset "Nodes" begin
    @testset "Empty Nodes" begin
        node = PhyNode()
        @test blisknown(node) == false
        @test confisknown(node) == false
        @test isempty(node) == true
    end
    @testset "Filled nodes (basic fields)" begin
        # create an empty node
        node = PhyNode("test", 10.0, 0.99)
        @test name(node) == "test"
        @test blisknown(node) == true
        @test branchlength(node) == 10.0
        @test confisknown(node) == true
        @test confidence(node) == 0.99
        @test isempty(node) == false
        # remove the fields
        confidence!(node, -1.0)
        @test confisknown(node) == false
        branchlength!(node, -1.0)
        @test blisknown(node) == false
        name!(node, "")
        @test isempty(node) == true
    end
    @testset "Relationships" begin
        @testset "Equality" begin
            a = PhyNode()
            b = PhyNode()
            @test isequal(a, b) == true
            name!(a, "a")
            @test isequal(a, b) == false
            name!(b, "a")
            @test isequal(a, b) == true
        end
        @testset "Parent-child" begin
            # we start with no relationship
            a = PhyNode()
            b = PhyNode()
            @test haschildren(a) == false
            @test hasparent(b) == false
            @test parentisself(b) == true
            @test isunlinked(a) == true
            @test countchildren(a) == 0
            @test haschild(a, b) == false
            # now create the parent-child relationship
            graft!(a, b)
            @test haschildren(a) == true
            @test hasparent(b) == true
            @test parentisself(b) == false
            @test isunlinked(a) == false
            @test countchildren(a) == 1
            @test haschild(a, b) == true
        end
        @testset "Root" begin
            a = PhyNode()
            @test isroot(a) == false
            b = PhyNode()
            graft!(a, b)
            @test isroot(a) == true
            c = PhyNode()
            graft!(c, a)
            @test isroot(a) == false
            @test isroot(c) == true
        end
        @testset "Siblings" begin
            a = PhyNode()
            b = PhyNode()
            graft!(a, b)
            @test length(siblings(b)) == 0
            c = PhyNode()
            graft!(a, c)
            @test length(siblings(b)) == 1
            @test in(c, siblings(b)) == true
        end
        @testset "Ancestor-descendant" begin
            a = PhyNode()
            b = PhyNode()
            @test isancestral(a, [b]) == false
            graft!(a, b)
            @test isancestral(a, [b]) == true
            c = PhyNode()
            graft!(a, c)
            @test isancestral(a, [c]) == true
        end
    end
end

## Relationships
#- build and destroy relationships between nodes

## Measures
#- distances, lengths, etc.


end # TestPhylo
