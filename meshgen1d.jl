
"Represents a 1D mesh with nodes and elements"
struct Mesh1D
    nnodes::Int   # Number of nodes
    node_coordinates::Array{Float64, 1}  # Coordinate of each node
    nelements::Int  # Number of elements

    # Map element ID to nodes
    element_nodes::Array{Int, 2}

    "Create a 1D mesh.
        node_coordinates is a 1D array of node locations
        elements is a 2D array, mapping elements to node indices"
    function Mesh1D(node_coordinates, elements)
        # Check that elements are in range 
        @assert maximum(elements) == length(node_coordinates)
        @assert minimum(elements) == 1
        
        new(length(node_coordinates),
            node_coordinates,
            size(elements)[1],
            elements)
    end
end

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
    
    return Mesh1D(VX, EToV)
end

