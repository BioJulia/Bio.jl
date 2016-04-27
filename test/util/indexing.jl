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
        newnames = [:Sixth, :Seventh, :Eighth, :Ninth, :Tenth]
        newnamestext = ["Sixth", "Seventh", "Eighth", "Ninth", "Tenth"]
        @test names!(i, newnames) == SingleIndex(newnames)
        @test_throws names!(i, newnames[1:4])
        @test names!(i, newnamestext) == SingleIndex(newnamestext)
        @test_throws names!(i, newnamestext[1:4])

    end

end




end
