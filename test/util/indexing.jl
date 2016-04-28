module TestIndexing

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Indexing

@testset "SingleIndex" begin
    symbolnames = [:First, :Second, :Third, :Fourth, :Fifth]
    textnames = ["First", "Second", "Third", "Fourth", "Fifth"]
    i = SingleIndex(symbolnames)
    e = SingleIndex()
    t = SingleIndex(textnames)
    @testset "Construction" begin
        @test e == SingleIndex()
        @test t != SingleIndex()
        @test i != SingleIndex()
        @test i == SingleIndex(convert(Vector{Symbol}, symbolnames))
        @test t == i
    end
    @testset "Basic Operators" begin
        @test length(i) == 5
        @test length(e) == 0
        @test length(t) == 5
        @test names(i) == symbolnames
        @test names(t) == symbolnames
        @test names(e) == Symbol[]
        @test isequal(copy(i), i) && !is(copy(i), i)
        @test isequal(copy(e), e) && !is(copy(e), e)
        @test isequal(copy(t), t) && !is(copy(t), t)
        for n in symbolnames
            @test haskey(i, n) == true
            @test haskey(t, n) == true
            @test haskey(e, n) == false
        end
        for n in textnames
            @test haskey(i, n) == true
            @test haskey(t, n) == true
            @test haskey(e, n) == false
        end
        for n in 1:5
            @test haskey(i, n) == true
            @test haskey(t, n) == true
            @test haskey(e, n) == false
        end
        @test keys(i) == symbolnames
        @test keys(t) == symbolnames
        @test keys(e) != symbolnames
        @test keys(e) == Symbol[]
    end
    @testset "Manipulation methods" begin
        @testset "Setting and resetting names" begin

            # Ok let's give i a load of new names
            newnames = [:Sixth, :Seventh, :Eighth, :Ninth, :Tenth]
            newnamestext = ["Sixth", "Seventh", "Eighth", "Ninth", "Tenth"]
            @test names!(i, newnames) == SingleIndex(newnames)
            @test_throws ArgumentError names!(i, newnames[1:4])
            @test names!(i, newnamestext) == SingleIndex(newnamestext)
            @test_throws ArgumentError names!(i, newnamestext[1:4])

            # Now let's manipulate some of those names by renaming them

            # Make a new index by replacing name Sixth with First,
            # and test that only :Sixth is different
            rn = rename(i, :Sixth, :First)
            @test (rn != SingleIndex(symbolnames)) && (rn != SingleIndex(newnames))
            @test rn.names[1] == :First
            @test rn.names[2] == :Seventh
            @test rn.names[3] == :Eighth
            @test rn.names[4] == :Ninth
            @test rn.names[5] == :Tenth

            # Let's make a new index by replacing the names, Sixth Eights, and
            # Tenth, with First, Third, and Fifth.
            rntwo = rename(i, [:Sixth, :Eighth, :Tenth], [:First, :Third, :Fifth])
            @test (rntwo != SingleIndex(symbolnames)) && (rntwo != SingleIndex(newnames))
            @test rn != rntwo

            # Now let's test that renaming worked correctly, some names should
            # match names in symbolnames, and some should match names in
            # newnames.
            @test (rntwo.names[1] != newnames[1]) && (rntwo.names[1] == symbolnames[1])
            @test (rntwo.names[2] == newnames[2]) && (rntwo.names[2] != symbolnames[2])
            @test (rntwo.names[3] != newnames[3]) && (rntwo.names[3] == symbolnames[3])
            @test (rntwo.names[4] == newnames[4]) && (rntwo.names[4] != symbolnames[4])
            @test (rntwo.names[5] != newnames[5]) && (rntwo.names[5] == symbolnames[5])

            @test (rntwo.names[1] != i.names[1]) && (rntwo.names[1] == t.names[1])
            @test (rntwo.names[2] == i.names[2]) && (rntwo.names[2] != t.names[2])
            @test (rntwo.names[3] != i.names[3]) && (rntwo.names[3] == t.names[3])
            @test (rntwo.names[4] == i.names[4]) && (rntwo.names[4] != t.names[4])
            @test (rntwo.names[5] != i.names[5]) && (rntwo.names[5] == t.names[5])
        end
        @testset "getindex" begin
            @test i[1] == 1
            @test i[[true, false, false, true, true]] == [1, 4, 5]
            @test i[1:3] == [1:3;]
            @test (i[:Sixth] == 1) && (i[:Tenth] == 5) && (i[:Eighth] == 3)
            @test (i["Sixth"] == 1) && (i["Tenth"] == 5) && (i["Eighth"] == 3)
            @test i[[1,3,5]] == [1,3,5]
            @test i[[:Seventh, :Eighth, :Ninth]] == [2, 3, 4]
            @test i[["Seventh", "Eighth", "Ninth"]] == [2, 3, 4]
        end
    end
end


end
