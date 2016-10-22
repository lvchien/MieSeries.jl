using MieSeries
using Base.Test

include("test_harmonics.jl")

# write your own tests here
from_matlab = [
    8.098707839783351
    2.064213995903995
    1.251251398359278
    1.045470773461690
    0.958170336867781
    0.905146609154803
    0.854560436057232
    0.818937250719909
    0.862677950996668
    0.834386518680957]


radius = 1.2
eps0 = 1.5
mu0 = 1.0
sphere = MieSeries.PECSphere(radius, eps0, mu0)
expnsn = MieSeries.excite(sphere, 2π, :left)

Θ = linspace(0.0,   π, 10)
Φ = linspace(0.0, 0.0,  1)
pts = [ [cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)] for θ in Θ for ϕ in Φ]

N = length(pts)
E = zeros(N)
H = zeros(N)
for (i,p) in enumerate(pts)
    Ex, Ey, Ez, Hx, Hy, Hz = MieSeries.farfield(sphere, expnsn, p)
    E[i] = norm([Ex,Ey,Ez])
    H[i] = norm([Ex,Ey,Ez])
end

for i in eachindex(E)
    @eval @test_approx_eq_eps E[$i] from_matlab[$i] sqrt(eps(Float64))
end

scattered_from_matlab = [
    0.334493148260097   0.448055290925219   0.474627521463338   0.321583691519947   0.217397836178558
    0.944332487407548                   0   0.982002821143344   0.397084866195122   0.239462161862970
    2.446818207854618                   0                   0   0.455389391557446   0.259887640573314
    0.944332487407548                   0   0.982002821143344   0.397084866195122   0.239462161862970
    0.334493148260097   0.448055290925219   0.474627521463338   0.321583691519947   0.217397836178558]

X = linspace(-2.0, 2.0, 5)
Z = linspace(-2.0, 4.0, 5)
P = [ [x, 0, z] for z in Z for x in X ]

A = zeros(Complex128, length(P))
for (i,p) in enumerate(P)
    norm(p) < 1.2 && continue
    ex, ey, ez, _ = scatteredfield(sphere, expnsn, p)
    A[i] = norm([ex,ey,ez])
end
A = reshape(A, length(X), length(Z))

for i in eachindex(A)
    @eval @test_approx_eq A[$i] scattered_from_matlab[$i]
end
