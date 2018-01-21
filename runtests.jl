#!/usr/bin/env julia

using Base.Test


tic()

for file in ["test-utils.jl", "test-grid1d.jl"]
    println("Tests in ", file)
    @time include(file)
end

toc()
