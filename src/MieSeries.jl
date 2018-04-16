module MieSeries

export PECSphere
export excite
export farfield
export scatteredfield

using SpecialFunctions

"""
    Ylm(K, r)

Computes the spherical harmonics up to order `K-1` in direction `r` .
"""
function Ylm(K, r)

    r = r / norm(r)
    x = r[1]; y = r[2]; z = r[3]

    T = eltype(r)
    Y = zeros(Complex{T}, K*K)

    Y[1] = 1/√4π

    facphip =  x + im*y;
    facphim = -x + im*y;
    for m ∈ 1:(K-1)
        Y[(m+1)*(m+1)] = -Y[m*m] * facphip * sqrt((m+0.5)/m)
        Y[m*m+1] = -Y[(m-1)*(m-1)+1] * facphim * sqrt((m+0.5)/m)
    end

    Y[3] = √3 * z * Y[1]
    for l ∈ 2:(K-1)
        Y[(l+1)*(l+1)-1] = sqrt(2l+1) * z * Y[l*l]
        Y[l*l+2] = sqrt(2*l+1) * z * Y[(l-1)*(l-1)+1]
    end

    for l ∈ 2:(K-1)
        for m ∈ (2-l):(l-2)
            Y[l*(l+1)+m+1] =
                sqrt((4*l^2-1)/(l^2-m^2)) * z * Y[l*(l-1)+m+1] -
                sqrt(((2*l+1)*((l-1)^2-m^2))/((l^2-m^2)*(2*l-3))) * Y[(l-1)*(l-2)+m+1]
        end
    end

    Y
end


function Xlm(K, r)

    r = r / norm(r)
    x = r[1]; y = r[2]; z = r[3];

    Y = Ylm(K, r);
    Xx = zeros(Y); Xy = zeros(Y); Xz = zeros(Y);
    Yx = zeros(Y); Yy = zeros(Y); Yz = zeros(Y);
    Zx = zeros(Y); Zy = zeros(Y); Zz = zeros(Y);

    for l = 1:(K-1)
        Ymin = sqrt(2*l) * Y[(l+1)*(l+1)-1];
        Xx[l*(l+1) + l + 1] = 1/sqrt(4*l*(l+1)) * Ymin;
        Xy[l*(l+1) + l + 1] = im/sqrt(4*l*(l+1)) * Ymin;
        Xz[l*(l+1) + l + 1] = l * Y[(l+1)*(l+1)] / sqrt(l*(l+1));

        # m = -l
        Yplus = sqrt(2*l) * Y[l*l + 2];
        Xx[l*(l+1) - l + 1] = 1/sqrt(4*l*(l+1)) * Yplus;
        Xy[l*(l+1) - l + 1] =-im/sqrt(4*l*(l+1)) * Yplus;
        Xz[l*(l+1) - l + 1] = -l * Y[l*l + 1] / sqrt(l*(l+1));

        # other values of m
        for m = (1-l):(l-1)
            Yplus = sqrt((l-m)*(l+m+1)) * Y[l*(l+1) + m+1 + 1];
            Ymin  = sqrt((l+m)*(l-m+1)) * Y[l*(l+1) + m-1 + 1];
            Xx[l*(l+1) + m + 1] = 1/sqrt(4*l*(l+1)) * (Yplus + Ymin);
            Xy[l*(l+1) + m + 1] =-im/sqrt(4*l*(l+1)) * (Yplus - Ymin);
            Xz[l*(l+1) + m + 1] = m * Y[l*(l+1) + m + 1] / sqrt(l*(l+1));
        end
    end

    # Y
    Yx = Y * x;
    Yy = Y * y;
    Yz = Y * z;

    # Z aka n × X
    Zx = y * Xz - z * Xy
    Zy = z * Xx - x * Xz
    Zz = x * Xy - y * Xx


    Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz
end


"""
    MNlm(numterms, position, wavenumber, besselfn)

Compute the multipole solutions to Maxwell's equations up to the orbital angular
momentum `numterms`. The argument `besselfn` can be one of `sbessel` or `shankel`.

The convention is such that the temporal variation of a fasor is

```math
a(t) = Re(a*e^{-iωt})
```

In particular, the Hankel function of the first kind corresponds to outgoing waves.

The multipole solutions are defined by

```math
    M_{l,m}(p) = besself(kr) X_{lm}(Θ,ϕ)
    N_{l,m}(p) = \frac{1}{k} ∇ × (besself(kr) X_{lm}(Θ,ϕ))
```

with ``X_{l,m}`` the vector spherical harmonics as defined in [Jackson, Classical Electrodynamics],
i.e. normalised in ``L^2(S)"`` with ``S`` the unit sphere.
"""
function MNlm(K, r, k, f)

    # Compute the spherical harmonics
    Xx, Xy, Xz, Yx, Yy, Yz, Px, Py, Pz = Xlm(K, r)

    # Argument van de besselfuncties (ook fysisch houden) en golfgetal:
    nor = k*norm(r)
    r = k * r

    # De functies invullen
    C = complex(eltype(r))
    Mx = zeros(C, K*K)
    My = zeros(C, K*K)
    Mz = zeros(C, K*K)
    Nx = zeros(C, K*K)
    Ny = zeros(C, K*K)
    Nz = zeros(C, K*K)

    for l = 1:(K-1)
        bessel = f(l, nor);
        besseldiff = ((l+1)*f(l-1, nor) - l*f(l+1, nor)) / (2*l+1);

        for m = -l:l
            lm = l*(l+1)+m+1
        	Mx[lm] = bessel * Xx[lm]
            My[lm] = bessel * Xy[lm]
            Mz[lm] = bessel * Xz[lm]

            factor = im * sqrt(l*(l+1)) * bessel / nor
            Nx[lm] = factor * Yx[lm] + besseldiff * Px[lm]
            Ny[lm] = factor * Yy[lm] + besseldiff * Py[lm]
            Nz[lm] = factor * Yz[lm] + besseldiff * Pz[lm]
        end
    end
    Mx, My, Mz, Nx, Ny, Nz
end


type PECSphere
    radius
    background_permittivity
    background_permeability
end

type Expansion
    numterms
    omega
    la1
    la2
    ra1
    ra2
end

"""
    excite(object, omega, pol)

Compute the expansion coefficients of the field traces in terms of harmonics.
The polarisation can be choosen from `:left`, `:right`, or ':x' (default).
"""
function excite(obj, omega, pol)

    r = obj.radius
    T = typeof(r)
    C = complex(T)
    K = 20

    # polarisation
    if pol == :left
        pols = [-1]
    elseif pol == :right
        pols = [+1]
    elseif pol == :x
        pols = [-1,1]
    end

    ra1 = zeros(C,K)
    ra2 = zeros(C,K)
    la1 = zeros(C,K)
    la2 = zeros(C,K)

    eps0 = obj.background_permittivity
    mu0  = obj.background_permeability
    k0 = omega * sqrt(eps0 * mu0);

    # Bessels, hankels and derivatives
    ll = collect(1:K)
    jk0   = sbessel.(ll, k0*r)
    djk0 = ((ll+1).*sbessel.(ll-1, k0*r)  - ll.*sbessel.(ll+1, k0*r))  ./ (2*ll+1)
    hk0   = shankel.(ll, k0*r)
    dhk0 = ((ll+1).*shankel.(ll-1, k0*r)  - ll.*shankel.(ll+1, k0*r))  ./ (2*ll+1)

    # compute the expansion coefficients for all requested polarisations
    for p in pols
        for l in ll
            a1 = (1*(-im)^l * sqrt(4*pi*(2*l+1)) * jk0[l]) / (-hk0[l])
            a2 = (p*(-im)^l * sqrt(4*pi*(2*l+1)) * djk0[l])/(-dhk0[l])
            # rescale for linear polarisation
            if pol == :x
                a1 = a1/2
                a2 = a2/2
            end
            # assign
            if p == -1
                la1[l] = a1
                la2[l] = a2
            elseif p == +1
                ra1[l] = a1
                ra2[l] = a2
            end
        end
    end

    Expansion(K, omega, la1, la2, ra1, ra2)
end


function farfield(obj, expansion, punten)

    # allocate space
    T = eltype(punten)
    C = complex(T)

    Ex = zero(C)
    Ey = zero(C)
    Ez = zero(C)

    Hx = zero(C)
    Hy = zero(C)
    Hz = zero(C)

    # basisfuncties in posities punten
    Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz = Xlm(expansion.numterms+1, punten)

    eps0 = obj.background_permittivity
    mu0  = obj.background_permeability
    k = expansion.omega * sqrt(eps0 * mu0);
    #k = expansion.omega * sqrt(obj.eps0*obj.mu0)
    Z = sqrt(mu0 / eps0)

    # build fields
    for p = [-1, 1]
        if p == -1
            a1 = expansion.la1
            a2 = expansion.la2
        else
            a1 = expansion.ra1
            a2 = expansion.ra2
        end
        m = -p
        for l in 1:expansion.numterms

            h =  im^(l+1) / k # shankel(l, k*r);
            dh = h * (-im)    # ((l+1).*shankel(l-1, k*r)  - l.*shankel(l+1, k*r))  ./ (2*l+1);

            # contributie in X_{l,m}
            Ex += a1[l] * h .* Xx[l*(l+1) + m + 1]
            Ey += a1[l] * h .* Xy[l*(l+1) + m + 1]
            Ez += a1[l] * h .* Xz[l*(l+1) + m + 1]
            Hx += im / Z * a2[l] * h .* Xx[l*(l+1) + m + 1]
            Hy += im / Z * a2[l] * h .* Xy[l*(l+1) + m + 1]
            Hz += im / Z * a2[l] * h .* Xz[l*(l+1) + m + 1]

            # contributie in Z_{l,m}
            Ex += a2[l] * dh .* Zx[l*(l+1) + m + 1]
            Ey += a2[l] * dh .* Zy[l*(l+1) + m + 1]
            Ez += a2[l] * dh .* Zz[l*(l+1) + m + 1]
            Hx += im/Z * a1[l] * dh .* Zx[l*(l+1) + m + 1]
            Hy += im/Z * a1[l] * dh .* Zy[l*(l+1) + m + 1]
            Hz += im/Z * a1[l] * dh .* Zz[l*(l+1) + m + 1]

        end
    end

    Ex, Ey, Ez, Hx, Hy, Hz
end



function scatteredfield(obj, expansion, punten)

    # outputarrays
    T = eltype(punten)
    C = complex(T)

    Ex = zero(C)
    Ey = zero(C)
    Ez = zero(C)

    Hx = zero(C)
    Hy = zero(C)
    Hz = zero(C)

    eps0 = obj.background_permittivity
    mu0  = obj.background_permeability
    k = expansion.omega * sqrt(eps0*mu0)
    Z = sqrt(mu0 / eps0)

    # multipole solutions of Maxwell's equations
    Mx, My, Mz, Nx, Ny, Nz = MNlm(expansion.numterms+1, punten, k, shankel)

    # radius
    r = norm(punten)

    # build fields
    for p in [-1, 1]
        if p == -1
            A1 = expansion.la1
            A2 = expansion.la2
        else
            A1 = expansion.ra1
            A2 = expansion.ra2
        end
        m = -p
        for l in 1:expansion.numterms

            # linear index
            lm = l*(l+1)+m+1

            h  = shankel(l, k*r)
            dh = ((l+1)*shankel(l-1, k*r)  - l*shankel(l+1, k*r))  / (2*l+1)

            a1 = A1[l]
            a2 = A2[l]

            # contribution in M_{l,m}
            Ex += a1 * Mx[lm]
            Ey += a1 * My[lm]
            Ez += a1 * Mz[lm]
            Hx += a2 * im / Z * Mx[lm]
            Hy += a2 * im / Z * My[lm]
            Hz += a2 * im / Z * Mz[lm]

            # contributie in N_{l,m}
            Ex += a2 * Nx[lm]
            Ey += a2 * Ny[lm]
            Ez += a2 * Nz[lm]
            Hx += a1 * im / Z * Mx[lm]
            Hy += a1 * im / Z * My[lm]
            Hz += a1 * im / Z * Mz[lm]

        end
    end

    Ex, Ey, Ez, Hx, Hy, Hz
end


"""
    sbessel(nu, x)

returns the spherical Bessel functions jnu(x).
"""
sbessel(nu, x) = sqrt(pi ./(2x)) .* besselj(nu + 0.5, x)

"""
    shankel(nu, x)

Computes the spherical hankel function (of the first kind).
"""
shankel(nu, x) = sqrt(pi ./(2x)) .* besselh(nu + 0.5, x)


end # module
