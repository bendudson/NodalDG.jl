# Routines for 2D grids

function Vandermonde(basis, quadrature)
    V = zeros(length(quadrature.locations), basis.order)
    for j=1:basis.order
        V[:,j] = basisFunction(basis, quadrature.locations, j)
    end
end

"Compute scaled warp function at order N based on rout
interpolation nodes"
function Warpfactor(basis, quadrature, rout)
    N = len(quadrature.locations)

    # Equidistant distribution
    eq = equidistant(N)

    # Compute V based on equidistant nodes
    Veq = Vandermonde(basis, eq)

    # Evaluate basis function at rout
    Pmat = zeros(basis.order, length(rout))
    for i=1:basis.order
        Pmat[i,:] = basisFunction(basis, rout, i)
    end
    Lmat = Veq' / Pmat

    # Compute warp factor
    warp = Lmat' * (quadrature.locations - req)

    # Scale factor
    
end

"Compute (x,y) nodes in equilateral triangle for polynomial
of order N"
function Nodes2D(N::Int)
    # Optimum correction to warp factor, to minimise
    # Lebesque constant for interpolant
    alpopt = [0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
              0.9800, 1.0999, 1.2832, 1.3648, 1.4773,
              1.4959, 1.5743, 1.5770, 1.6223, 0.6258]

    # Set optimised parameter, alpha, depending on order
    alpha = N < 16 ? alpopt[N] : 5.0 / 3

    # Total number of nodes
    Np = (N + 1) * (N + 2) / 2

    # Create equidistributed nodes on equilateral triangle
    L1 = zeros(Np, 1)
    L3 = zeros(Np, 1)
    sk = 1
    for n = 1:N+1
        for m = 1:N+2-n
            L1[sk] = (n - 1) / N
            L3[sk] = (m - 1) / N
            sk += 1
        end
    end
    L2 = 1.0 - L1 - L3
    x = -L2 + L3
    y = (-L2 - L3 + 2 * L1)/sqrt(3.0)

    # Compute blending function at each node for each edge
    blend1 = 4 * L2 * L3
    blend2 = 4 * L1 * L3
    blend3 = 4 * L1 * L2

    # Amount of warp for each node, for each edge
    warpf1 = Warpfactor(N, L3 - L2)
    warpf2 = Warpfactor(N, L1 - L3)
    warpf3 = Warpfactor(N, L2 - L1)

    
end
