
module Testing

import ..NodalDG

using Test

# Simple test of JacobiP
@test NodalDG.JacobiP(0.0, 0, 0, 0) == [1.0 /sqrt(2.0)]

# Test Gauss-Lobatto points

r, w = NodalDG.gausslobatto(3)
@test isapprox(r, [-1, 0, 1])
@test isapprox(w, [1/3, 4/3, 1/3])

r,w = NodalDG.gausslobatto(4)
@test isapprox(r, [-1, -sqrt(1/5), sqrt(1/5), 1])
@test isapprox(w, [1/6, 5/6, 5/6, 1/6])

# Tests based on p48 in textbook, Fig 3.1
# This tests JacobiP and LGL points

r, w = NodalDG.gausslobatto(7)

# Fig 3.1d
Vref = [0.707107 -1.22474 1.58114 -1.87083 2.12132 -2.34521 2.54951;
        0.707107 -1.01681 0.844182 -0.346644 -0.278373 0.807539 -1.05741;
        0.707107 -0.57422 -0.269222 0.833675 -0.504704 -0.365166 0.846707;
        0.707107 0.0 -0.790569 -0.0 0.795495 0.0 -0.796722;
        0.707107 0.57422 -0.269222 -0.833675 -0.504704 0.365166 0.846707;
        0.707107 1.01681 0.844182 0.346644 -0.278373 -0.807539 -1.05741;
        0.707107 1.22474 1.58114 1.87083 2.12132 2.34521 2.54951]

@test isapprox( NodalDG.Vandermonde1D(6,r), Vref, atol=1e-5 )

# Tests based on p54, Fig 3.4

# order N = 1
r, w = NodalDG.gausslobatto(2)
V = NodalDG.Vandermonde1D(1, r)

Dref1 = [-0.5 0.5;
         -0.5 0.5]
@test isapprox( NodalDG.Dmatrix1D(1, r, V), Dref1 )

# order N = 2
r, w = NodalDG.gausslobatto(3)
V = NodalDG.Vandermonde1D(2, r)

Dref2 = [-1.5 2.0 -0.5;
         -0.5 0.0 0.5;
         0.5 -2.0 1.5]

@test isapprox( NodalDG.Dmatrix1D(2, r, V), Dref2 )

# order N = 4
r, w = NodalDG.gausslobatto(5)
V = NodalDG.Vandermonde1D(4, r)

Dref4 = [-5.0 6.7565 -2.66667 1.41016 -0.5;
         -1.24099 2.51054e-16 1.74574 -0.763763 0.25901;
         0.375 -1.33658 -0.0 1.33658 -0.375;
         -0.25901 0.763763 -1.74574 -6.56798e-16 1.24099;
         0.5 -1.41016 2.66667 -6.7565 5.0]

@test maximum(abs, NodalDG.Dmatrix1D(4, r, V) - Dref4) < 1e-5

end
