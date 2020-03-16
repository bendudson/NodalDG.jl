
"Represents a 1D mesh with vertices and elements"
struct Mesh1D
    nvertices::Int   # Number of vertices
    vertex_coordinates::Array{Float64, 1}  # Coordinate of each vertex
    nelements::Int  # Number of elements

    # Map element ID to vertices
    element_vertices::Array{Int, 2}

    "Create a 1D mesh.
        vertex_coordinates is a 1D array of node locations
        elements is a 2D array, mapping elements to vertex indices"
    function Mesh1D(vertex_coordinates, elements)
        # Check that elements are in range 
        @assert maximum(elements) == length(vertex_coordinates)
        @assert minimum(elements) == 1
        
        new(length(vertex_coordinates),
            vertex_coordinates,
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

