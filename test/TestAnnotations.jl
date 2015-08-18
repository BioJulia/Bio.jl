module TestAnnotations

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Annotations

firstField() = collect(1.0:5.0)
secondField() = [false, true, false, true, true]
thirdField() = collect(1:5)
fourthField() = ["Pugh", "Barney McGrew", "Cuthbert", "Dibble", "Grub"]

firstConstructorTest() = AnnotationContainer(tuple(firstField(), secondField()), Tuple{:a, :b})
secondConstructorTest() = AnnotationContainer(tuple(firstField(), secondField(), fourthField()), Tuple{:a, :b, :c})

@testset "Annotations") begin
    @testset "Construction") begin
        @test typeof(firstConstructorTest()) == AnnotationContainer{Tuple{Array{Float64, 1}, Array{Bool, 1}}, Tuple{:a, :b}}
        @test @annots(:a = firstField(), :b = secondField()) == firstConstructorTest()
        @test typeof(secondConstructorTest()) == AnnotationContainer{Tuple{Array{Float64, 1}, Array{Bool, 1}, Array{ASCIIString, 1}}, Tuple{:a, :b, :c}}
        @test @annots(:a = firstField(), :b = secondField(), :c = fourthField()) == secondConstructorTest()
        @test AnnotationContainer(firstConstructorTest(), thirdField(), Field{:c}) == @annots(:a = firstField(), :b = secondField(), :c = thirdField())
    end

    @testset "Fetching Values") begin
        @test typeof(firstConstructorTest[Field{:a}]) == Array{Float64, 1}
        @test firstConstructorTest[Field{:a}] == firstField
        @test typeof(firstConstructorTest[Field{:b}]) == Array{Bool, 1}
        @test firstConstructorTest[Field{:b}] == secondField
        @test typeof(secondConstructorTest[Field{:a}]) == Array{Float64, 1}
        @test secondConstructorTest[Field{:a}] == firstField
        @test typeof(secondConstructorTest[Field{:b}]) == Array{Bool, 1}
        @test secondConstructorTest[Field{:b}] == secondField
        @test typeof(secondConstructorTest[Field{:c}]) == Array{ASCIIString, 1}
        @test secondConstructorTest[Field{:c}] == fourthField

        @test typeof(firstConstructorTest[1, Field{:a}]) == Float64
        @test typeof(firstConstructorTest[2, Field{:a}]) == Float64
        @test typeof(firstConstructorTest[3, Field{:a}]) == Float64
        @test typeof(firstConstructorTest[4, Field{:a}]) == Float64
        @test typeof(firstConstructorTest[5, Field{:a}]) == Float64
        @test firstConstructorTest[1, Field{:a}] == 1.0
        @test firstConstructorTest[2, Field{:a}] == 2.0
        @test firstConstructorTest[3, Field{:a}] == 3.0
        @test firstConstructorTest[4, Field{:a}] == 4.0
        @test firstConstructorTest[5, Field{:a}] == 5.0
        @test typeof(firstConstructorTest[1, Field{:b}]) == Bool
        @test typeof(firstConstructorTest[2, Field{:b}]) == Bool
        @test typeof(firstConstructorTest[3, Field{:b}]) == Bool
        @test typeof(firstConstructorTest[4, Field{:b}]) == Bool
        @test typeof(firstConstructorTest[5, Field{:b}]) == Bool
        @test firstConstructorTest[1, Field{:b}] == false
        @test firstConstructorTest[2, Field{:b}] == true
        @test firstConstructorTest[3, Field{:b}] == false
        @test firstConstructorTest[4, Field{:b}] == true
        @test firstConstructorTest[5, Field{:b}] == true

        @test typeof(secondConstructorTest[1, Field{:a}]) == Float64
        @test typeof(secondConstructorTest[2, Field{:a}]) == Float64
        @test typeof(secondConstructorTest[3, Field{:a}]) == Float64
        @test typeof(secondConstructorTest[4, Field{:a}]) == Float64
        @test typeof(secondConstructorTest[5, Field{:a}]) == Float64
        @test secondConstructorTest[1, Field{:a}] == 1.0
        @test secondConstructorTest[2, Field{:a}] == 2.0
        @test secondConstructorTest[3, Field{:a}] == 3.0
        @test secondConstructorTest[4, Field{:a}] == 4.0
        @test secondConstructorTest[5, Field{:a}] == 5.0
        @test typeof(secondConstructorTest[1, Field{:b}]) == Bool
        @test typeof(secondConstructorTest[2, Field{:b}]) == Bool
        @test typeof(secondConstructorTest[3, Field{:b}]) == Bool
        @test typeof(secondConstructorTest[4, Field{:b}]) == Bool
        @test typeof(secondConstructorTest[5, Field{:b}]) == Bool
        @test secondConstructorTest[1, Field{:b}] == false
        @test secondConstructorTest[2, Field{:b}] == true
        @test secondConstructorTest[3, Field{:b}] == false
        @test secondConstructorTest[4, Field{:b}] == true
        @test secondConstructorTest[5, Field{:b}] == true
        @test typeof(secondConstructorTest[1, Field{:c}]) == ASCIIString
        @test typeof(secondConstructorTest[2, Field{:c}]) == ASCIIString
        @test typeof(secondConstructorTest[3, Field{:c}]) == ASCIIString
        @test typeof(secondConstructorTest[4, Field{:c}]) == ASCIIString
        @test typeof(secondConstructorTest[5, Field{:c}]) == ASCIIString
        @test secondConstructorTest[1, Field{:c}] == "Pugh"
        @test secondConstructorTest[2, Field{:c}] == "Barney McGrew"
        @test secondConstructorTest[3, Field{:c}] == "Cuthbert"
        @test secondConstructorTest[4, Field{:c}] == "Dibble"
        @test secondConstructorTest[5, Field{:c}] == "Grub"
    end
end

end
