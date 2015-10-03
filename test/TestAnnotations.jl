module TestAnnotations

using FactCheck, Bio.Annotations

const firstField = collect(1.0:5.0)
const secondField = [false, true, false, true, true]
const thirdField = collect(1:5)
const fourthField = ["Pugh", "Barney McGrew", "Cuthbert", "Dibble", "Grub"]

facts("Annotations") do
    context("Construction") do
        firstConstructorTest = AnnotationContainer(tuple(firstField, secondField), Tuple{:a, :b})
        secondConstructorTest = AnnotationContainer(tuple(firstField, secondField, fourthField), Tuple{:a, :b, :c})
        @fact typeof(firstConstructorTest) --> AnnotationContainer{Tuple{Array{Float64, 1}, Array{Bool, 1}}, Tuple{:a, :b}}
        @fact @annotations(:a = firstField, :b = secondField) --> firstConstructorTest
        @fact typeof(secondConstructorTest) --> AnnotationContainer{Tuple{Array{Float64, 1}, Array{Bool, 1}, Array{ASCIIString, 1}}, Tuple{:a, :b, :c}}
        @fact @annotations(:a = firstField, :b = secondField, :c = fourthField) --> secondConstructorTest
        @fact AnnotationContainer(firstConstructorTest, thirdField, Field{:c}) --> @annotations(:a = firstField, :b = secondField, :c = thirdField)
    end

    context("Fetching Values") do
        @fact typeof(firstConstructorTest[Field{:a}]) --> Array{Float64, 1}
        @fact firstConstructorTest[Field{:a}] --> firstField
        @fact typeof(firstConstructorTest[Field{:b}]) --> Array{Bool, 1}
        @fact firstConstructorTest[Field{:b}] --> secondField
        @fact typeof(secondConstructorTest[Field{:a}]) --> Array{Float64, 1}
        @fact secondConstructorTest[Field{:a}] --> firstField
        @fact typeof(secondConstructorTest[Field{:b}]) --> Array{Bool, 1}
        @fact secondConstructorTest[Field{:b}] --> secondField
        @fact typeof(secondConstructorTest[Field{:c}]) --> Array{ASCIIString, 1}
        @fact secondConstructorTest[Field{:c}] --> fourthField

        @fact typeof(firstConstructorTest[1, Field{:a}]) --> Float64
        @fact typeof(firstConstructorTest[2, Field{:a}]) --> Float64
        @fact typeof(firstConstructorTest[3, Field{:a}]) --> Float64
        @fact typeof(firstConstructorTest[4, Field{:a}]) --> Float64
        @fact typeof(firstConstructorTest[5, Field{:a}]) --> Float64
        @fact firstConstructorTest[1, Field{:a}] --> 1.0
        @fact firstConstructorTest[2, Field{:a}] --> 2.0
        @fact firstConstructorTest[3, Field{:a}] --> 3.0
        @fact firstConstructorTest[4, Field{:a}] --> 4.0
        @fact firstConstructorTest[5, Field{:a}] --> 5.0
        @fact typeof(firstConstructorTest[1, Field{:b}]) --> Bool
        @fact typeof(firstConstructorTest[2, Field{:b}]) --> Bool
        @fact typeof(firstConstructorTest[3, Field{:b}]) --> Bool
        @fact typeof(firstConstructorTest[4, Field{:b}]) --> Bool
        @fact typeof(firstConstructorTest[5, Field{:b}]) --> Bool
        @fact firstConstructorTest[1, Field{:b}] --> false
        @fact firstConstructorTest[2, Field{:b}] --> true
        @fact firstConstructorTest[3, Field{:b}] --> false
        @fact firstConstructorTest[4, Field{:b}] --> true
        @fact firstConstructorTest[5, Field{:b}] --> true

        @fact typeof(secondConstructorTest[1, Field{:a}]) --> Float64
        @fact typeof(secondConstructorTest[2, Field{:a}]) --> Float64
        @fact typeof(secondConstructorTest[3, Field{:a}]) --> Float64
        @fact typeof(secondConstructorTest[4, Field{:a}]) --> Float64
        @fact typeof(secondConstructorTest[5, Field{:a}]) --> Float64
        @fact secondConstructorTest[1, Field{:a}] --> 1.0
        @fact secondConstructorTest[2, Field{:a}] --> 2.0
        @fact secondConstructorTest[3, Field{:a}] --> 3.0
        @fact secondConstructorTest[4, Field{:a}] --> 4.0
        @fact secondConstructorTest[5, Field{:a}] --> 5.0
        @fact typeof(secondConstructorTest[1, Field{:b}]) --> Bool
        @fact typeof(secondConstructorTest[2, Field{:b}]) --> Bool
        @fact typeof(secondConstructorTest[3, Field{:b}]) --> Bool
        @fact typeof(secondConstructorTest[4, Field{:b}]) --> Bool
        @fact typeof(secondConstructorTest[5, Field{:b}]) --> Bool
        @fact secondConstructorTest[1, Field{:b}] --> false
        @fact secondConstructorTest[2, Field{:b}] --> true
        @fact secondConstructorTest[3, Field{:b}] --> false
        @fact secondConstructorTest[4, Field{:b}] --> true
        @fact secondConstructorTest[5, Field{:b}] --> true
        @fact typeof(secondConstructorTest[1, Field{:c}]) --> ASCIIString
        @fact typeof(secondConstructorTest[2, Field{:c}]) --> ASCIIString
        @fact typeof(secondConstructorTest[3, Field{:c}]) --> ASCIIString
        @fact typeof(secondConstructorTest[4, Field{:c}]) --> ASCIIString
        @fact typeof(secondConstructorTest[5, Field{:c}]) --> ASCIIString
        @fact secondConstructorTest[1, Field{:c}] --> "Pugh"
        @fact secondConstructorTest[2, Field{:c}] --> "Barney McGrew"
        @fact secondConstructorTest[3, Field{:c}] --> "Cuthbert"
        @fact secondConstructorTest[4, Field{:c}] --> "Dibble"
        @fact secondConstructorTest[5, Field{:c}] --> "Grub"
    end
end

end
