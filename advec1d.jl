#
# Advection in 1D
#

using NodalDG

# Low storage Runge-Kutta coefficients
rk4a = [ 0.0 ...
         -567301805773.0/1357537059087.0 ...
         -2404267990393.0/2016746695238.0 ...
         -3550918686646.0/2091501179385.0  ...
         -1275806237668.0/842570457699.0]
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0]
rk4c = [ 0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0]

"Evaluate RHS flux in 1D advection"
function AdvecRHS1D(grid, u, time, a)

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
    
    # Runge-Kutta residual storage  
    resu = zeros(grid.Np, grid.K); 

    # compute time step size
    xmin = minimum(abs, grid.x[1,:]-grid.x[2,:])
    CFL=0.75
    dt = CFL/a*xmin
    dt = .5*dt
    Nsteps = ceil(FinalTime/dt)
    dt = FinalTime/Nsteps
    
    # outer time step loop 
    for tstep=1:Nsteps
        resu[:] = 0.0
        for INTRK = 1:5
            timelocal = time + rk4c[INTRK]*dt
            rhsu = AdvecRHS1D(grid, u, timelocal, a)
            resu = rk4a[INTRK]*resu + dt*rhsu
            u += rk4b[INTRK]*resu
        end
        # Increment time
        time += dt
    end
    return u
end

# Order of polymomials used for approximation 
N = 3
K = 12

# Generate simple mesh
Nv, VX, K, EToV = NodalDG.MeshGen1D(0.0,2.0*pi,K)

# Initialize solver and construct grid and metric
grid = NodalDG.Grid1D(VX, EToV, N)

using PyPlot


# Set initial conditions
u = sin.(grid.x)

# Solve Problem
outputDt = 0.1
time = 0.0
for i = 1:10
    println("Time = ", time)
    plot(grid.x, u)
    u = Advec1D(grid, u, time, time + outputDt)
    time += outputDt
end

show()
