module Testing

import ..NodalDG

using Test

@test_throws AssertionError NodalDG.Mesh1D([0, 1, 2],
                                           [[1 2]; [2 3]; [3 4]])

@test_throws AssertionError NodalDG.Mesh1D([0, 1, 2],
                                           [[0 1]; [1 2]])

let mesh = NodalDG.Mesh1D([0, 1, 2, 3],
                          [[1 2]; [2 3]; [3 4]])
    @test mesh.nvertices == 4
    @test isapprox(mesh.vertex_coordinates, [0, 1, 2, 3])
    @test mesh.nelements == 3
    @test mesh.element_vertices == [[1 2]; [2 3]; [3 4]]
end

end
