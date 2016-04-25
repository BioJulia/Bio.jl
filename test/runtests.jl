# all test targets
available_targets = [
    "align",
    "phylo",
    "intervals",
    "seq",
    "services",
    "structure",
    "tools"
]

if isempty(ARGS)
    # run all available test targets
    targets = available_targets
else
    targets = ARGS
    invalids = setdiff(targets, available_targets)
    if !isempty(invalids)
        error("there are invalid test targets: ", join(invalids, ", "))
    end
end

include("TestFunctions.jl")

for target in targets
    include("$target/runtests.jl")
end

for target in readdir("util")
    include("$target")
end
