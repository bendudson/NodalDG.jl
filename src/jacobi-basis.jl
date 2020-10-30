
using SpecialFunctions: gamma

"
Evaluate Jacobi polynomial of type 
(alpha, beta) > -1, (alpha+beta != -1)
at points x for order N. Returns array of values, same length as x

See: https://en.wikipedia.org/wiki/Jacobi_polynomials

NOTE: This calculates the normalised Jacobi polynomial

"
function JacobiP(x, alpha, beta, N::Int)

    nvalues = length(x)
    
    PL = zeros(Float64, nvalues, N+1)

    gamma0 = 2^(alpha+beta+1)/(alpha+beta+1) * gamma(alpha+1) * gamma(beta+1)/gamma(alpha+beta+1)
    PL[:,1] .= 1.0 / sqrt(gamma0)
    
    if N == 0
        return PL[:,1]
    end

    gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0
    PL[:,2] = ((alpha+beta+2)*x/2 .+ (alpha-beta)/2)/sqrt(gamma1)

    if N == 1
        return PL[:,2]
    end

    aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3))

    # Forward recurrence
    for i = 1:(N-1)
        h1 = 2*i + alpha + beta
        anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3))
        bnew = - (alpha^2 - beta^2)/h1/(h1+2)
        PL[:,i+2] = (1/anew)*(-aold * PL[:,i] + (x .- bnew) .* PL[:,i+1])
        aold = anew
    end
    return PL[:,end]
end

"Evaluate the derivative of the Jacobi polynomial
of type (alpha,beta) > -1 at points r for order N

## Returns
array of the same length as r
"
function GradJacobiP(r, alpha, beta, N::Int)
    if N == 0
        zeros(length(r))
    else
        sqrt(N * (N+alpha+beta+1))*JacobiP(r, alpha+1, beta+1, N-1)
    end
end

struct JacobiBasis1D
    order::Int  # Number of polynomials
end

function basisFunction(basis::JacobiBasis1D, locations, index::Int)
    JacobiP(locations, 0, 0, index-1)
end

function basisGradient(basis::JacobiBasis1D, locations, index::Int)
    GradJacobiP(locations,0,0,index-1)
end
