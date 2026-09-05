module MidpointExtrapolation
using LinearAlgebra
export skewextract, midpoint_K, cayley, invcayley, extrapolated_step
skewextract(B)=tril(B,-1)-transpose(tril(B,-1))
function cayley(K)
    n=size(K,1); I0=Matrix{eltype(K)}(I,n,n); (I0-K)\(I0+K)
end
function invcayley(C)
    n=size(C,1); I0=Matrix{eltype(C)}(I,n,n); (C-I0)/(C+I0)
end
function midpoint_K(G,delta; K0=nothing, tol=1e-12, maxiter=20)
    n=size(G,1); I0=Matrix{eltype(G)}(I,n,n)
    K = K0===nothing ? (delta/2)*skewextract(G) : copy(K0)
    for r in 1:maxiter
        B=(I0+K) \ (G/(I0-K))
        Kn=(delta/2)*skewextract(B)
        if norm(Kn-K,Inf) <= tol*(1+norm(Kn,Inf)); return Kn,r; end
        K=Kn
    end
    return K,maxiter
end
"""One macrostep. Aof(t,Qbase) should return A(t) in physical coordinates."""
function extrapolated_step(Qn,t,h,Aof; tol=1e-12,maxiter=20)
    # half 1
    G1=transpose(Qn)*Aof(t+h/4,Qn)*Qn
    Kh1,r21=midpoint_K(G1,h/2;tol=tol,maxiter=maxiter); C1=cayley(Kh1); Qh=Qn*C1
    # half 2
    G2=transpose(Qh)*Aof(t+3h/4,Qh)*Qh
    Kh2,r22=midpoint_K(G2,h/2;tol=tol,maxiter=maxiter); C2=cayley(Kh2)
    Cnet=C1*C2; K2=invcayley(Cnet)
    # full, predicted by K2
    Gf=transpose(Qn)*Aof(t+h/2,Qn)*Qn
    K1,r1=midpoint_K(Gf,h;K0=K2,tol=tol,maxiter=maxiter)
    Kext=(4K2-K1)/3; Qnew=Qn*cayley(Kext)
    err=norm((K2-K1)/3, Frobenius())
    return Qnew,(r1=r1,r21=r21,r22=r22,error=err,K1=K1,K2=K2,Kext=Kext)
end
end
