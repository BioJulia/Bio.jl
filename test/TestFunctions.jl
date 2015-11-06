module TestFunctions


export get_bio_fmt_specimens,
    random_array


function get_bio_fmt_specimens()
    path = Pkg.dir("Bio", "test", "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone --depth 1 https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
end

function random_array(n::Integer, elements, probs)
    cumprobs = cumsum(probs)
    x = Array(eltype(elements), n)
    for i in 1:n
        x[i] = elements[searchsorted(cumprobs, rand()).start]
    end
    return x
end


end
