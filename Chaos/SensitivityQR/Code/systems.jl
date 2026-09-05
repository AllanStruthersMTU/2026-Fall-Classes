module SensitivitySystems
using LinearAlgebra
export rossler!, lorenz!, hyperrossler!, lorenz96!, rossler_jac!, lorenz_jac!, hyperrossler_jac!, lorenz96_jac!

function rossler!(du,u,p,t)
    a,b,c = p
    x,y,z = u
    du[1] = -y-z; du[2] = x+a*y; du[3] = b+z*(x-c)
end
function rossler_jac!(A,u,p,t)
    a,b,c = p; x,y,z=u
    A .= 0; A[1,2]=-1; A[1,3]=-1; A[2,1]=1; A[2,2]=a; A[3,1]=z; A[3,3]=x-c
end
function lorenz!(du,u,p,t)
    sigma,rho,beta=p; x,y,z=u
    du[1]=sigma*(y-x); du[2]=x*(rho-z)-y; du[3]=x*y-beta*z
end
function lorenz_jac!(A,u,p,t)
    sigma,rho,beta=p; x,y,z=u
    A .= 0; A[1,1]=-sigma; A[1,2]=sigma; A[2,1]=rho-z; A[2,2]=-1; A[2,3]=-x; A[3,1]=y; A[3,2]=x; A[3,3]=-beta
end
function hyperrossler!(du,u,p,t)
    a,b,c,d=p; x,y,z,w=u
    du[1]=-y-z; du[2]=x+a*y+w; du[3]=b+x*z; du[4]=-c*z+d*w
end
function hyperrossler_jac!(A,u,p,t)
    a,b,c,d=p; x,y,z,w=u
    A .= 0; A[1,2]=-1; A[1,3]=-1; A[2,1]=1; A[2,2]=a; A[2,4]=1; A[3,1]=z; A[3,3]=x; A[4,3]=-c; A[4,4]=d
end
function lorenz96!(du,u,F,t)
    N=length(u)
    @inbounds for i in 1:N
        im2=mod1(i-2,N); im1=mod1(i-1,N); ip1=mod1(i+1,N)
        du[i]=(u[ip1]-u[im2])*u[im1]-u[i]+F
    end
end
function lorenz96_jac!(A,u,F,t)
    N=length(u); A .= 0
    @inbounds for i in 1:N
        im2=mod1(i-2,N); im1=mod1(i-1,N); ip1=mod1(i+1,N)
        A[i,ip1]+=u[im1]; A[i,im2]+=-u[im1]; A[i,im1]+=u[ip1]-u[im2]; A[i,i]+=-1
    end
end
end
