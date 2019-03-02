#
# Advection in 1D
#

include("NodalDG.jl")

using .NodalDG

"Evaluate RHS flux in 1D advection"
function AdvecRHS1D(u, time, grid, a)

    # form field differences at faces
    alpha = 0.5
    
    du = zeros(grid.Nfp*grid.Nfaces,grid.K); 
    du[:] = (u[grid.vmapM[:]]-u[grid.vmapP[:]]).*(a*grid.nx[:]-(1-alpha)*abs.(a*grid.nx[:]))/2
    
    # impose boundary condition at x=0
    uin = -sin(a*time)
    du[grid.mapI] = (u[grid.vmapI]- uin ).*(a*grid.nx[grid.mapI]-(1-alpha)*abs.(a*grid.nx[grid.mapI]))/2
    du[grid.mapO] = 0

    # compute right hand sides of the semi-discrete PDE
    rhsu = -a*grid.rx.*(grid.Dr*u) + grid.LIFT*(grid.Fscale.*(du))
    
    return rhsu
end

"Integrate 1D advection until FinalTime starting with
initial the condition, u"
function Advec1D(grid, u, time, FinalTime)
    # advection speed
    a = 2*pi 

    # compute time step size
    xmin = minimum(abs, grid.x[1,:]-grid.x[2,:])
    CFL=0.75
    dt = CFL/a*xmin
    dt = .5*dt

    return NodalDG.lserk!(AdvecRHS1D, u,
                          time, FinalTime, dt,
                          args=(grid, a))
end

# Order of polymomials used for approximation 
N = 5
K = 5

# Generate simple mesh
mesh = NodalDG.MeshGen1D(0.0,2.0*pi,K)

# Initialize solver and construct grid and metric
grid = NodalDG.Grid1D(mesh, N)

using PyPlot: plot, show


# Set initial conditions
state = sin.(grid.x)

# Solve Problem
outputDt = 0.1
sim_time = 0.0

for i = 1:10
    global state, sim_time
    println("Time = ", sim_time)
    plot(grid.x, state)
    state = Advec1D(grid, state, sim_time, sim_time + outputDt)
    sim_time += outputDt
end

show()
