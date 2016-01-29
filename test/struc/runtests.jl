module TestStruc

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Struc

include("test_model.jl")
include("test_pdb.jl")
include("test_spatial.jl")

end
