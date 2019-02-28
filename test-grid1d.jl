module Testing

using ..NodalDG
using Test

VX = [0.0, 1.0, 2.0, 3.0]
EToV = [1 2; 2 3; 3 4]

g = NodalDG.Grid1D(VX, EToV, 2)

println(g.vmapM, g.vmapP)


# Test based on p58/9 in book. Note errata on website
EToV = [1 2; 2 3; 3 5; 5 4]


end

