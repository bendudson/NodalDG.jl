
using LinearAlgebra
using SparseArrays


"Construct the Element-To-Element (EToE) and
Element-To-Face (EtoF) maps
"
function Connect1D(EToV)
    K = size(EToV)[1] # Number of cells   
    Nfaces = 2
    
    TotalFaces = Nfaces*K  # Total number of faces
    Nv = K+1  # Number of vertices
    
    # Local face to local vertex connections
    vn = [1,2]
    
    FToVrows = 1:TotalFaces
    FToVcols = zeros(TotalFaces)
    FToVvals = ones(Int, TotalFaces)
    
    row = 1
    for k=1:K
        for face=1:Nfaces
            FToVcols[row] = EToV[k,vn[face]]
            row += 1
        end
    end
    FToV = sparse(FToVrows, FToVcols, FToVvals)
    
    # Global face to global face sparse array
    FToF = FToV * FToV' - I
    
    faces1, faces2, _ = findnz(FToF) # Find non-zero elements
    
    
    # Convert face global number to element and face numbers
    element1 = floor.(Int, (faces1 .- 1)/Nfaces) .+ 1
    face1 = mod.( faces1 .- 1, Nfaces ) .+ 1
    element2 = floor.(Int, (faces2 .- 1)/Nfaces) .+ 1
    face2 = mod.( faces2 .- 1, Nfaces ) .+ 1

    # Rearrange into Nelements x Nfaces sized arrays

    # This initialises all faces to connect to themselves
    # Faces which connect to other cells are then overwritten
    EToE = (1:K) * ones(Int, 1, Nfaces)
    EToF = ones(Int, K,1)*(1:Nfaces)'

    # this index code doesn't work as expected
    #ind = sub2ind([K,Nfaces], element1, face1)
    #EToE[ind] = element2
    #EToF[ind] = face2

    for i = 1:length(element1)
        EToE[element1[i], face1[i]] = element2[i]
        EToF[element1[i], face1[i]] = face2[i]
    end
    
    return EToE, EToF
end

"
K    Number of cells
Nfp  Number of nodes per face
Np   Number of nodes per cell
"
function BuildMaps1D(K::Int, Nfaces::Int, Nfp::Int, Np::Int, Fmask, EToE, EToF, x)
    
    nodeids = reshape(1:K*Np, Np, K)
    
    vmapM   = zeros(Int, Nfp, Nfaces, K) 
    vmapP   = zeros(Int, Nfp, Nfaces, K)

    for k1=1:K
        for f1=1:Nfaces
            # find index of face nodes with respect to volume node ordering
            vmapM[:,f1,k1] .= nodeids[Fmask[f1], k1]
        end
    end

    for k1=1:K
        for f1=1:Nfaces
            # find neighbor
            k2 = EToE[k1,f1]
            f2 = EToF[k1,f1]

            # find volume node numbers of left and right nodes 
            vidM = vmapM[:,f1,k1]
            vidP = vmapM[:,f2,k2]
            
            x1  = x[vidM]
            x2  = x[vidP]
    
            # Compute distance matrix
            D = sum(abs2, x1 -x2)
            if D < 1e-5
                vmapP[:,f1,k1] = vidP
            end
        end
    end
    
    #vmapP = vmapP[:]
    #vmapM = vmapM[:]

    # Find where vmapP == vmapM
    mapB = findall(iszero, vmapP - vmapM)
    vmapB = vmapM[mapB]

    # Create specific left (inflow) and right (outflow) maps
    mapI = 1
    mapO = K*Nfaces
    vmapI = 1
    vmapO = K*Np;
    
    return vmapP, vmapM, mapB, vmapB, mapI, mapO, vmapI, vmapO
end

# Compute LIFT
function Lift1D(Np, Nfaces, Nfp, V)
    Emat = zeros(Np, Nfaces*Nfp)
    
    # Define Emat
    Emat[1,1] = 1.0
    Emat[Np,2] = 1.0
    
    # inv(mass matrix)*\s_n (L_i,L_j)_{edge_n}
    V*(V'*Emat);
end


"Describes a 1D grid with K cells"
mutable struct Grid1D
    K::Int # Number of cells
    Np::Int # Number of points per cell
    Nfp::Int # Number of points on a face
    Nfaces::Int # Number of faces per cell

    VX::Array{Float64,1} # Array of vertex locations (K+1)
    EToV::Array{Int,2}   # Cell vertices for each cell

    r::Array{Float64,1} # node locations in reference cell [-1,1]
    w::Array{Float64,1}

    V    # Vandermonde matrix on reference cell
    invV # Inverse Vandermonde matrix
    Dr   # Differentiation matrix

    LIFT # Surface integral terms
    
    J    # Jacobian [point, cell]
    rx   # Inverse Jacobian
    
    x::Array{Float64, 2} # Physical coordinates of the grid points [point,cell]

    Fscale
    
    nx  # Outward pointing normal at element faces

    vmapP # Map of external points for each cell face
    vmapM # Map of internal points for each cell face
    
    mapB
    vmapB

    mapI
    mapO
    vmapI
    vmapO
    
    function Grid1D(vertices::Array{Float64,1}, cells::Array{Int,2}, N::Int)

        # Create an instance. Not yet initialised
        # This is to avoid a long constructor parameter list
        g = new()
        
        g.VX = vertices
        g.EToV = cells

        g.K = size(g.EToV)[1] # Number of cells
        g.Np = N+1 # Number of points per cell (including edges)
        g.Nfaces = 2 # Number of faces per cell
        g.Nfp = 1  # Number of points per face
        
        # Compute LGL grid on reference cell
        g.r, g.w = gausslobatto(N+1)

        # Build reference element matrices
        g.V = Vandermonde1D(N, g.r)
        g.invV = inv(g.V)

        # Differentiation matrix on reference element
        g.Dr = Dmatrix1D(N, g.r, g.V)

        g.LIFT = Lift1D(g.Np, g.Nfaces, g.Nfp, g.V)
        
        # Coordinates of all the nodes
        va = g.EToV[:,1]'
        vb = g.EToV[:,2]'
        g.x = ones(N+1,1)*g.VX[va] .+ 0.5*(g.r .+ 1)*(g.VX[vb] - g.VX[va])

        # Geometric factors
        g.J = g.Dr * g.x
        g.rx = 1.0 ./ g.J

        fmask1 = 1
        fmask2 = N+1
        Fmask = [fmask1;fmask2]
        Fx = g.x[Fmask,:]  # Locations of cell edges for each cell

        g.Fscale = 1.0 ./ g.J[Fmask,:]
        
        # Compute outward pointing normals
        g.nx = zeros(g.Nfaces, g.K)
        g.nx[1,:] .= -1
        g.nx[2,:] .= 1

        # Build global connectivity arrays
        EToE, EToF = Connect1D(g.EToV)
        
        g.vmapP, g.vmapM, g.mapB, g.vmapB, g.mapI, g.mapO, g.vmapI, g.vmapO = BuildMaps1D(g.K, 2, 1, g.Np, Fmask, EToE, EToF, g.x)

        return g
    end
end
