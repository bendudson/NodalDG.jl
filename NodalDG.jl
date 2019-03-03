
# Top-level source for NodalDG
#
#
# Dependencies: FastGaussQuadrature  ( https://github.com/ajt60gaibb/FastGaussQuadrature.jl )
#


module NodalDG

include("meshgen1d.jl")

include("grid1d.jl")

include("lserk.jl")

include("jacobi-basis.jl")

end

