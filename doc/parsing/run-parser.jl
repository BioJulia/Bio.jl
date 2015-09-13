
include("bed3.jl")

for entry in open("data.bed", BED3)
    @show entry
end

