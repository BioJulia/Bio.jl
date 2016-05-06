module TestIndexing

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Indexing

@testset "Indexer" begin
    symbolnames = [:First, :Second, :Third, :Fourth, :Fifth]
    textnames = ["First", "Second", "Third", "Fourth", "Fifth"]
    groups = UnitRange{UInt}[1:5, 5:10, 10:15, 15:20, 20:25]
    groups_eight = UnitRange{UInt8}[1:5, 5:10, 10:15, 15:20, 20:25]
    int_index_one = Indexer(symbolnames)
    int_index_two = Indexer(textnames, UInt64)
    int_index_three = Indexer(symbolnames, UInt8)
    g_index_one = Indexer(symbolnames, groups)
    g_index_two = Indexer(textnames, groups)
    g_index_three = Indexer(textnames, groups_eight)
    empty_index = Indexer(UInt64)

    @testset "Correct Construction" begin
        @test int_index_one == int_index_two
        @test int_index_one == int_index_three
        @test int_index_two == int_index_three
        @test g_index_one == g_index_two
        @test typeof(int_index_one) == typeof(int_index_two)
        @test typeof(int_index_one) != typeof(int_index_three)
        @test typeof(int_index_two) != typeof(int_index_three)
        @test typeof(g_index_one) == typeof(g_index_two)
        @test typeof(g_index_one) != typeof(g_index_three)
        @test typeof(g_index_two) != typeof(g_index_three)
        @test typeof(int_index_one) == Indexer{UInt64}
        @test typeof(int_index_two) == Indexer{UInt64}
        @test typeof(int_index_three) == Indexer{UInt8}
        @test typeof(g_index_one) == Indexer{UnitRange{UInt64}}
        @test typeof(g_index_two) == Indexer{UnitRange{UInt64}}
        @test typeof(g_index_three) == Indexer{UnitRange{UInt8}}
        @test typeof(empty_index) == Indexer{UInt64}
    end

    @testset "Basic Operators" begin
        @test length(int_index_one) == 5
        @test length(int_index_two) == 5
        @test length(int_index_three) == 5
        @test length(g_index_one) == 5
        @test length(g_index_two) == 5
        @test length(g_index_three) == 5
        @test names(int_index_one) == symbolnames
        @test names(int_index_two) == symbolnames
        @test names(int_index_three) == symbolnames
        @test names(g_index_one) == symbolnames
        @test names(g_index_two) == symbolnames
        @test names(g_index_three) == symbolnames
        @test names(empty_index) == Symbol[]
        @test isequal(copy(int_index_one), int_index_one) && !is(copy(int_index_one), int_index_one)
        @test isequal(copy(int_index_two), int_index_two) && !is(copy(int_index_two), int_index_two)
        @test isequal(copy(int_index_three), int_index_three) && !is(copy(int_index_three), int_index_three)
        @test isequal(copy(g_index_one), g_index_one) && !is(copy(g_index_one), g_index_one)
        @test isequal(copy(g_index_two), g_index_two) && !is(copy(g_index_two), g_index_two)
        @test isequal(copy(g_index_three), g_index_three) && !is(copy(g_index_three), g_index_three)
        @test isequal(copy(empty_index), empty_index) && !is(copy(empty_index), empty_index)
        for n in symbolnames
            @test haskey(int_index_one, n) == true
            @test haskey(int_index_two, n) == true
            @test haskey(int_index_three, n) == true
            @test haskey(g_index_one, n) == true
            @test haskey(g_index_two, n) == true
            @test haskey(g_index_three, n) == true
            @test haskey(empty_index, n) == false
        end
        for n in textnames
            @test haskey(int_index_one, n) == true
            @test haskey(int_index_two, n) == true
            @test haskey(int_index_three, n) == true
            @test haskey(g_index_one, n) == true
            @test haskey(g_index_two, n) == true
            @test haskey(g_index_three, n) == true
            @test haskey(empty_index, n) == false
        end
        for n in 1:5
            @test haskey(int_index_one, n) == true
            @test haskey(int_index_two, n) == true
            @test haskey(int_index_three, n) == true
            @test haskey(g_index_one, n) == true
            @test haskey(g_index_two, n) == true
            @test haskey(g_index_three, n) == true
            @test haskey(empty_index, n) == false
        end
        @test keys(int_index_one) == symbolnames
        @test keys(int_index_two) == symbolnames
        @test keys(int_index_three) == symbolnames
        @test keys(g_index_one) == symbolnames
        @test keys(g_index_two) == symbolnames
        @test keys(g_index_three) == symbolnames
        @test keys(empty_index) == Symbol[]
        @test keys(empty_index) != symbolnames
    end

    @testset "Manipulation methods" begin
        @testset "Setting and resetting names" begin

            # Ok let's give int_index_one a load of new names
            newnames = [:Sixth, :Seventh, :Eighth, :Ninth, :Tenth]
            newnamestext = ["Sixth", "Seventh", "Eighth", "Ninth", "Tenth"]
            @test names!(int_index_one, newnames) == Indexer(newnames)
            @test_throws ArgumentError names!(int_index_one, newnames[1:4])
            @test names!(int_index_two, newnamestext) == Indexer(newnamestext)
            @test_throws ArgumentError names!(int_index_two, newnamestext[1:4])

            # Now let's manipulate some of those names by renaming them

            # Make a new index by replacing name :Sixth with :First,
            # and test that only :Sixth is different
            rn = rename(int_index_one, :Sixth, :First)
            @test (rn != Indexer(symbolnames)) && (rn != Indexer(newnames))
            @test names(rn)[1] == :First
            @test names(rn)[2] == :Seventh
            @test names(rn)[3] == :Eighth
            @test names(rn)[4] == :Ninth
            @test names(rn)[5] == :Tenth

            # Let's modify int_index_one by replacing the names, Sixth Eights, and
            # Tenth, with First, Third, and Fifth.
            rename!(int_index_one, [:Sixth, :Eighth, :Tenth], [:First, :Third, :Fifth])
            @test int_index_one != Indexer(symbolnames)
            @test rntwo != Indexer(newnames)
            @test rn != int_index_one

            # Now let's test that renaming worked correctly, some names should
            # match names in symbolnames, and some should match names in
            # newnames.
            @test names(int_index_one)[1] != newnames[1]
            @test names(int_index_one)[1] == symbolnames[1]
            @test names(int_index_one)[2] == newnames[2]
            @test names(int_index_one)[2] != symbolnames[2]
            @test names(int_index_one)[3] != newnames[3]
            @test names(int_index_one)[3] == symbolnames[3]
            @test names(int_index_one)[4] == newnames[4]
            @test names(int_index_one)[4] != symbolnames[4]
            @test names(int_index_one)[5] != newnames[5]
            @test names(int_index_one)[5] == symbolnames[5]
        end
        @testset "getindex" begin
            @test int_index_one[1] == 1
            @test int_index_one[[true, false, false, true, true]] == [1, 4, 5]
            @test int_index_one[1:3] == [1:3;]
            @test int_index_one[:Sixth] == 1
            @test int_index_one[:Tenth] == 5
            @test int_index_one[:Eighth] == 3
            @test int_index_one["Sixth"] == 1
            @test int_index_one["Tenth"] == 5
            @test int_index_one["Eighth"] == 3
            @test int_index_one[[1,3,5]] == [1,3,5]
            @test int_index_one[[:Seventh, :Eighth, :Ninth]] == [2, 3, 4]
            @test int_index_one[["Seventh", "Eighth", "Ninth"]] == [2, 3, 4]
        end
    end
end

end
