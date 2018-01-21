
"Generate simple equidistant grid with K elements"
function MeshGen1D(xmin, xmax, K)

    # Number of vertices
    Nv = K+1
    
    # Generate node coordinates
    VX = zeros(Nv)
    for i = 1:Nv
        VX[i] = (xmax-xmin)*(i-1)/(Nv-1) + xmin
    end

    # read element to node connectivity
    EToV = zeros(Int, K, 2)
    for k = 1:K
        EToV[k,1] = k
        EToV[k,2] = k+1
    end
    
    return Nv, VX, K, EToV
end

