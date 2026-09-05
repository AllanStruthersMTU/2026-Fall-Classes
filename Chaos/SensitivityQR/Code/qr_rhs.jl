module ContinuousQR
using LinearAlgebra
export skewextract!, direct_q_rhs!
function skewextract!(Omega,B)
    n=size(B,1); fill!(Omega,0)
    @inbounds for j in 1:n, i in j+1:n
        Omega[i,j]=B[i,j]; Omega[j,i]=-B[i,j]
    end
    Omega
end
"""Direct continuous-Q RHS for a supplied A(t). State is vec(Q)."""
function direct_q_rhs!(dq,q,p,t)
    n,Afun=p; Q=reshape(q,n,n)
    A=Afun(t); H=A*Q; B=transpose(Q)*H; Omega=zeros(eltype(Q),n,n); skewextract!(Omega,B)
    dQ=Q*Omega; dq .= vec(dQ)
end
end
