
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
    print(Np,", ", Nfaces, ", ", Nfp)
    Emat = zeros(Np, Nfaces*Nfp)
    
    # Define Emat
    Emat[1,1] = 1.0
    Emat[Np,2] = 1.0
    
    # inv(mass matrix)*\s_n (L_i,L_j)_{edge_n}
    V*(V'*Emat);
end

struct Quadrature1D
    # locations in reference cell [-1,1]
    locations::Array{Float64, 1}
    weights::Array{Float64, 1}
end

"Compute LGL quadrature points on reference cell [-1,1]"
function gaussLobatto1D(npoints)
    r, w = gausslobatto(npoints)
    Quadrature1D(r, w)
end

"Describes a 1D grid with K cells"
mutable struct Grid1D
    
    mesh::Mesh1D # Defines node locations and cells
    
    npoints::Int # Number of points per cell, including edges (Np)
    npoints_per_face::Int # Number of points on a face
    nfaces::Int # Number of faces per cell

    quadrature::Quadrature1D
    
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
    
    function Grid1D(mesh::Mesh1D, quadrature::Quadrature1D)
        
        # Create an instance. Not yet initialised
        # This is to avoid a long constructor parameter list
        g = new()
        
        g.mesh = mesh
        
        K = mesh.nelements # Number of cells

        g.npoints = length(quadrature.locations)
        g.nfaces = 2 # Number of faces per cell
        g.npoints_per_face = 1  # Number of points per face

        g.quadrature = quadrature

        N = length(quadrature.locations)-1
        
        # Build reference element matrices
        g.V = Vandermonde1D(N, quadrature.locations)
        g.invV = inv(g.V)

        # Differentiation matrix on reference element
        g.Dr = Dmatrix1D(N, quadrature.locations, g.V)

        g.LIFT = Lift1D(g.npoints, g.nfaces, g.npoints_per_face, g.V)
        
        # Coordinates of all the nodes
        va = mesh.element_nodes[:,1]'
        vb = mesh.element_nodes[:,2]'
        g.x = ones(N+1,1)*mesh.node_coordinates[va] .+ 0.5*(quadrature.locations  .+ 1)*(mesh.node_coordinates[vb] - g.mesh.node_coordinates[va])

        # Geometric factors
        g.J = g.Dr * g.x
        g.rx = 1.0 ./ g.J

        fmask1 = 1
        fmask2 = N+1
        Fmask = [fmask1;fmask2]
        Fx = g.x[Fmask,:]  # Locations of cell edges for each cell

        g.Fscale = 1.0 ./ g.J[Fmask,:]
        
        # Compute outward pointing normals
        g.nx = zeros(g.nfaces, mesh.nelements)
        g.nx[1,:] .= -1
        g.nx[2,:] .= 1

        # Build global connectivity arrays
        EToE, EToF = Connect1D(mesh.element_nodes)
        
        g.vmapP, g.vmapM, g.mapB, g.vmapB, g.mapI, g.mapO, g.vmapI, g.vmapO = BuildMaps1D(mesh.nelements, 2, 1, g.npoints, Fmask, EToE, EToF, g.x)

        return g
    end
end
