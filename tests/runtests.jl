#!/usr/bin/env julia

using Test

using NodalDG.jl

for file in ["test-utils.jl", "test-mesh1d.jl", "test-grid1d.jl"]
    println("Tests in ", file)
    @time include(file)
end

