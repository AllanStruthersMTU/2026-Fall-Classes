module CayleyTools
using LinearAlgebra
export cayley, invcayley, local_cayley_rhs, skewextract
skewextract(B) = tril(B,-1)-transpose(tril(B,-1))
function cayley(K)
    n=size(K,1); I0=Matrix{eltype(K)}(I,n,n)
    return (I0-K) \ (I0+K)   # commuting polynomials in K; same Cayley map
end
function invcayley(C)
    n=size(C,1); I0=Matrix{eltype(C)}(I,n,n)
    return (C-I0)/(C+I0)
end
"""Step-local Cayley RHS K' for fixed Qn and a callback A(t)."""
function local_cayley_rhs(K,Qn,A)
    C=cayley(K); Q=Qn*C; B=transpose(Q)*(A*Q); Om=skewextract(B)
    return 0.5*((I+K)*Om*(I-K))
end
end
