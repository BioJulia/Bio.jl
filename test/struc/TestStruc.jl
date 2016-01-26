using Bio.Struc

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

include("test_model.jl")
include("test_pdb.jl")
include("test_spatial.jl")
