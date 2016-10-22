using MieSeries

a = 1.0
ϵ = 1.0
μ = 1.0
sphere = PECSphere(a, ϵ, μ)

ω   = 2π
k   = ω * √(ϵ * μ)
pol = :x
exc = excite(sphere, ω, pol)

Θ = linspace(0.0,   π, 10)
Φ = linspace(0.0, 0.0,  1)
P = [ [cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)] for θ in Θ for ϕ in Φ]

N = length(P)
E = zeros(N)
H = zeros(N)

for (i,p) in enumerate(P)
    Ex, Ey, Ez, Hx, Hy, Hz = farfield(sphere, exc, p)
    E[i] = norm([Ex,Ey,Ez])
    H[i] = norm([Ex,Ey,Ez])
end

# Define the grid on which to compute the scattered field
x = linspace(-2.0, 2.0, 3)
z = linspace(-4.0, 4.0, 400)
P = [ [ξ, 0, ζ] for ζ in z for ξ in x ]


# Compute the scattered field
Ex = zeros(Complex128, length(P))
Ey = zeros(Complex128, length(P))
Ez = zeros(Complex128, length(P))
for (i,p) in enumerate(P)
    norm(p) < a && continue
    ex, ey, ez, _ = scatteredfield(sphere, exc, p)
    Ex[i] = ex + exp(im*k*p[3])
    Ey[i] = ey
    Ez[i] = ex
end
Ex = reshape(Ex, length(x), length(z))
Ey = reshape(Ey, length(x), length(z))
Ez = reshape(Ez, length(x), length(z))

# Plot the field
A = real(sqrt(Ex.*conj(Ex)))

# using Plots
# n = div(length(x),2) + 1
# plot(A[n,:], m=:point)
# savefig("mie.svg")
