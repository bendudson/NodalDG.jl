
# Top-level source for NodalDG
#
#
# Dependencies: FastGaussQuadrature  ( https://github.com/ajt60gaibb/FastGaussQuadrature.jl )
#


module NodalDG
include("utils.jl")

include("grid1d.jl")

include("meshgen1d.jl")

include("lserk.jl")

end

