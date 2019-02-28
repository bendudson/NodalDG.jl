#!/usr/bin/env julia

using Test

include("NodalDG.jl")

for file in ["test-utils.jl", "test-grid1d.jl"]
    println("Tests in ", file)
    @time include(file)
end

