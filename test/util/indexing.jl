module TestIndexers

using Base.Test

using Bio.Indexers

@testset "Indexers" begin

    symbolnames = [:First, :Second, :Third, :Fourth, :Fifth]
    textnames = ["First", "Second", "Third", "Fourth", "Fifth"]
    groups = UnitRange{Int}[1:5, 6:10, 11:15, 16:20, 21:25]
    groups_eight = UnitRange{UInt8}[1:5, 6:10, 11:15, 16:20, 21:25]

    int_idxes = (Indexer(symbolnames),
                 Indexer(textnames, UInt64),
                 Indexer(symbolnames, UInt8))

    g_idxes = (Indexer(symbolnames, groups),
               Indexer(textnames, groups),
               Indexer(textnames, groups_eight))

    empty_index = Indexer(UInt64)

    @testset "Constructor Behaviour" for idxset in (int_idxes, g_idxes)
        for idxer in idxset
            @test idxer.names !== symbolnames && idxer.names == symbolnames
        end
    end

    @testset "Index equality" begin
        @testset "Single Indexes" for idxer in int_idxes
            for idx in int_idxes
                @test idxer == idx
            end
            for idx in g_idxes
                @test idxer != g_idxes
            end
            @test idxer != empty_index
        end
        @testset "Group Indexes" for idxer in g_idxes
            for idx in g_idxes
                @test idxer == idx
            end
            for idx in int_idxes
                @test idxer != idx
            end
            @test idxer != empty_index
        end
    end

    @testset "Basic Operators" begin
        # Tests for single indexes and group indexes.
        for idxset in (int_idxes, g_idxes)
            for idxer in idxset
                @test length(idxer) == 5
                @test names(idxer) == symbolnames
                @test names(idxer) !== symbolnames
                @test names(idxer) == idxer.names
                @test names(idxer) !== idxer.names
                for i in 1:5
                    @test haskey(idxer, symbolnames[i]) == true
                    @test haskey(idxer, textnames[i]) == true
                    @test haskey(idxer, i) == true
                end
                @test keys(idxer) == symbolnames
                @test keys(idxer) !== symbolnames
                idxcopy = copy(idxer)
                @test isequal(idxcopy, idxer)
                @test is(idxcopy, idxer) == false
                @test idxcopy.names == idxer.names
                @test idxcopy.names !== idxer.names
                @test idxcopy.lookup == idxer.lookup
                @test idxcopy.lookup !== idxer.lookup
            end
        end
        # Tests for empty index seperate because they are a bit different.
        @test names(empty_index) == Symbol[]
        @test names(empty_index) == Symbol[]
        @test names(empty_index) !== Symbol[]
        @test names(empty_index) !== empty_index.names
        @test names(empty_index) == empty_index.names
        for i in 1:5
            @test haskey(empty_index, symbolnames[i]) == false
            @test haskey(empty_index, textnames[i]) == false
            @test haskey(empty_index, i) == false
        end
        @test keys(empty_index) == Symbol[]
        @test keys(empty_index) !== Symbol[]
        idxcopy = copy(empty_index)
        @test isequal(idxcopy, empty_index)
        @test is(idxcopy, empty_index) == false
        @test idxcopy.names == empty_index.names
        @test idxcopy.names !== empty_index.names
        @test idxcopy.lookup == empty_index.lookup
        @test idxcopy.lookup !== empty_index.lookup
    end

    @testset "Manipulation methods" begin
        @testset "Setting and resetting names" begin

            # Ok let's give int_index_one a load of new names
            newnames = [:Sixth, :Seventh, :Eighth, :Ninth, :Tenth]
            newnamestext = ["Sixth", "Seventh", "Eighth", "Ninth", "Tenth"]
            @test names!(int_idxes[1], newnames) == Indexer(newnames)
            @test_throws ArgumentError names!(int_idxes[1], newnames[1:4])
            @test names!(int_idxes[2], newnamestext) == Indexer(newnamestext)
            @test_throws ArgumentError names!(int_idxes[2], newnamestext[1:4])

            # Now let's manipulate some of those names by renaming them

            # Make a new index by replacing name :Sixth with :First,
            # and test that only :Sixth is different
            rn = rename(int_idxes[1], :Sixth, :First)
            @test (rn != Indexer(symbolnames)) && (rn != Indexer(newnames))
            @test names(rn)[1] == :First
            @test names(rn)[2] == :Seventh
            @test names(rn)[3] == :Eighth
            @test names(rn)[4] == :Ninth
            @test names(rn)[5] == :Tenth

            # Let's modify int_index_one by replacing the names, Sixth Eights, and
            # Tenth, with First, Third, and Fifth.
            rename!(int_idxes[1], [:Sixth, :Eighth, :Tenth], [:First, :Third, :Fifth])
            @test int_idxes[1] != Indexer(symbolnames)
            @test rn != int_idxes[1]
            @test rn !== int_idxes[1]

            # Now let's test that renaming worked correctly, some names should
            # match names in symbolnames, and some should match names in
            # newnames.
            @test names(int_idxes[1])[1] != newnames[1]
            @test names(int_idxes[1])[1] == symbolnames[1]
            @test names(int_idxes[1])[2] == newnames[2]
            @test names(int_idxes[1])[2] != symbolnames[2]
            @test names(int_idxes[1])[3] != newnames[3]
            @test names(int_idxes[1])[3] == symbolnames[3]
            @test names(int_idxes[1])[4] == newnames[4]
            @test names(int_idxes[1])[4] != symbolnames[4]
            @test names(int_idxes[1])[5] != newnames[5]
            @test names(int_idxes[1])[5] == symbolnames[5]
        end
        @testset "getindex" begin
            idxer_one = int_idxes[1]
            @testset "single names" begin
                @test idxer_one[:First] == 1
                @test idxer_one[:Fifth] == 5
                @test idxer_one[:Third] == 3
                @test idxer_one["First"] == 1
                @test idxer_one["Fifth"] == 5
                @test idxer_one["Third"] == 3
                @test g_idxes[1][:First] == 1:5
                @test g_idxes[3][:First] == 1:5
            end
            @testset "vector of names" begin
                @test idxer_one[[:Seventh, :Third, :Ninth]] == [2, 3, 4]
                @test idxer_one[["Seventh", "Third", "Ninth"]] == [2, 3, 4]
                @test g_idxes[1][[:First, :Second]] == UnitRange{Int}[1:5, 6:10]
            end
            @testset "single integer" begin
                @test idxer_one[1] == 1
                @test g_idxes[1][4] == 16:20
            end
            @testset "multiple integers" begin
                @test idxer_one[1:3] == [1:3;]
                @test idxer_one[[1,3,5]] == [1,3,5]
            end
            @testset "booleans" begin
                @test idxer_one[[true, false, false, true, true]] == UInt[1, 4, 5]
            end
        end
    end
end

end
