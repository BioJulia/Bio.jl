module TestPhylo

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Phylo

@testset "Phylo" begin

    @testset "PhyNodes" begin
        @testset "Individual nodes" begin
            @testset "Constructors" begin
                @testset "Empty Nodes" begin
                    node = PhyNode()
                    @test has_branchlength(node) == false
                    @test has_support(node) == false
                    @test isempty(node) == true
                    @test node.name == ""
                    @test typeof(node) == PhyNode{Void, Void, Void}
                end
                @testset "Named nodes" begin
                    named = PhyNode("MRCA")
                    @test has_branchlength(named) == false
                    @test has_support(named) == false
                    @test isempty(named) == false
                    @test named.name == "MRCA"
                    @test typeof(named) == PhyNode{Void, Void, Void}
                end
                @testset "Nodes with support and branch" begin
                    nometa_one = PhyNode("MRCA1", 0.342, 0.982)
                    @test has_branchlength(nometa_one) == true
                    @test has_support(nometa_one) == true
                    @test isempty(nometa_one) == false
                    @test nometa_one.name == "MRCA1"
                    @test typeof(nometa_one) == PhyNode{Float64, Float64, Void}

                    nometa_two = PhyNode("MRCA2", true, 0.982)
                    @test_throws ErrorException has_branchlength(nometa_two)
                    @test has_support(nometa_two) == true
                    @test isempty(nometa_two) == false
                    @test nometa_two.name == "MRCA2"
                    @test typeof(nometa_two) == PhyNode{Bool, Float64, Void}

                    nometa_three = PhyNode("MRCA3", 'A', nothing)
                    @test_throws ErrorException has_branchlength(nometa_three)
                    @test has_support(nometa_three) == false
                    @test isempty(nometa_three) == false
                    @test nometa_three.name == "MRCA3"
                    @test typeof(nometa_three) == PhyNode{Char, Void, Void}
                end
            end
            @testset "Basic node fields" begin
                # Create a typical Phylogenetic Node with just branchlength
                # and bootstrap value.
                phylo_with_support = PhyNode("typical_one", 10.0, 0.99)
                @test name(phylo_with_support) == "typical_one"
                @test has_branchlength(phylo_with_support) == true
                @test branchlength(phylo_with_support) == 10.0
                @test has_support(phylo_with_support) == true
                @test support(phylo_with_support) == 0.99
                @test isempty(phylo_with_support) == false
                # remove the fields
                support!(phylo_with_support, -1.0)
                @test has_support(phylo_with_support) == false
                branchlength!(phylo_with_support, -1.0)
                @test has_branchlength(phylo_with_support) == false
                name!(phylo_with_support, "")
                @test isempty(phylo_with_support) == true
            end
            @testset "Node Relationships" begin
        #     @testset "Equality" begin
        #         a = PhyNode()
        #         b = PhyNode()
        #         @test isequal(a, b) == true
        #         name!(a, "a")
        #         @test isequal(a, b) == false
        #         name!(b, "a")
        #         @test isequal(a, b) == true
        #     end
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
                    # # now create the parent-child relationship
                    # graft!(a, b)
                    # @test haschildren(a) == true
                    # @test hasparent(b) == true
                    # @test parentisself(b) == false
                    # @test isunlinked(a) == false
                    # @test countchildren(a) == 1
                    # @test haschild(a, b) == true
                end
        #     @testset "Root" begin
        #         a = PhyNode()
        #         @test isroot(a) == false
        #         b = PhyNode()
        #         graft!(a, b)
        #         @test isroot(a) == true
        #         c = PhyNode()
        #         graft!(c, a)
        #         @test isroot(a) == false
        #         @test isroot(c) == true
        #     end
        #     @testset "Siblings" begin
        #         a = PhyNode()
        #         b = PhyNode()
        #         graft!(a, b)
        #         @test length(siblings(b)) == 0
        #         c = PhyNode()
        #         graft!(a, c)
        #         @test length(siblings(b)) == 1
        #         @test in(c, siblings(b)) == true
        #     end
        #     @testset "Ancestor-descendant" begin
        #         a = PhyNode()
        #         b = PhyNode()
        #         @test isancestral(a, [b]) == false
        #         graft!(a, b)
        #         @test isancestral(a, [b]) == true
        #         c = PhyNode()
        #         graft!(a, c)
        #         @test isancestral(a, [c]) == true
        #     end
            end
        end
    end
end




end # TestPhylo
