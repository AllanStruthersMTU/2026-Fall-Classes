using LinearAlgebra, DifferentialEquations, OrdinaryDiffEq, CSV, DataFrames, Plots
include("systems.jl"); include("qr_rhs.jl"); include("cayley_tools.jl"); include("midpoint_extrapolation.jl")
using .SensitivitySystems, .ContinuousQR, .CayleyTools, .MidpointExtrapolation

# Draft experiment driver. Production version should integrate y and Q consistently
# (or use a tight dense state solution y(t) as a common A(t) source for all Q methods).
function state_reference(f!,u0,p,tspan; abstol=1e-13,reltol=1e-13)
    prob=ODEProblem(f!,u0,tspan,p)
    solve(prob,Vern9();abstol=abstol,reltol=reltol,dense=true)
end
function make_Aof(sol,jac!,p,n)
    A=zeros(n,n)
    t -> begin
        jac!(A,sol(t),p,t); copy(A)
    end
end
function run_direct_Q(name,sol,jac!,p,n,tspan,alg; abstol=1e-8,reltol=1e-8)
    Aof=make_Aof(sol,jac!,p,n)
    q0=vec(Matrix{Float64}(I,n,n))
    rhs! = (dq,q,pp,t) -> ContinuousQR.direct_q_rhs!(dq,q,(n,Aof),t)
    prob=ODEProblem(rhs!,q0,tspan)
    @timed solve(prob,alg;abstol=abstol,reltol=reltol,save_everystep=true)
end
function demo()
    tspan=(0.0,100.0); u0=[1.0,0.0,0.0]; p=(0.2,0.2,5.7)
    sol=state_reference(rossler!,u0,p,tspan)
    oz=run_direct_Q("Rossler",sol,rossler_jac!,p,3,tspan,OwrenZen4())
    ts=run_direct_Q("Rossler",sol,rossler_jac!,p,3,tspan,Tsit5())
    println("OwrenZen4 direct-Q seconds = ",oz.time)
    println("Tsit5 direct-Q seconds = ",ts.time)
end
if abspath(PROGRAM_FILE)==@__FILE__; demo(); end
