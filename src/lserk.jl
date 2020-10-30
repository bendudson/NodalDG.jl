#

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


"
Integrates the RHS function rhsfunc
from t0 to t1, starting with solution u

This version modifies u 
"
function lserk!(rhsfunc, u, t0, t1, dt; args=())
    Nsteps = ceil((t1 - t0)/dt)
    dt = (t1 - t0)/Nsteps

    # Runge-Kutta residual storage  
    resu = similar(u)
    
    time = t0
    # outer time step loop 
    for tstep=1:Nsteps
        for INTRK = 1:5
            timelocal = time + rk4c[INTRK]*dt
            rhsu = rhsfunc(u, timelocal, args...)
            resu = rk4a[INTRK]*resu + dt*rhsu
            u += rk4b[INTRK]*resu
        end
        # Increment time
        time += dt
    end
    
    return u
end

"
Integrates the RHS function rhsfunc
from t0 to t1, starting with solution u0
"
function lserk(rhsfunc, u0, t0, t1, dt; args=())
    u = copy(u0)
    return lserk!(rhsfunc, u, t0, t1, dt; args=args)
end
