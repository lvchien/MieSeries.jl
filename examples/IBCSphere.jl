using MieSeries;

struct IBCSphere
	ZsOverEta0  # normalized impedance
    radius
    background_permittivity
    background_permeability
end

struct NdeltaSphere
	ZsOverEta0  # normalized impedance
	delta
    radius
    background_permittivity
    background_permeability
end

function relationSH(obj :: IBCSphere, ll)
	Zn = obj.ZsOverEta0

	vnxX = Zn*ones(ll);
	vX = Zn*ones(ll);
	return (vnxX, vX);
end
	
function relationSH(obj :: NdeltaSphere, ll)
	delta = obj.delta
	r = obj.radius
	Zn = obj.ZsOverEta0

	id = im/delta;
	jid   = MieSeries.sbessel(ll, id*r)     # Spherical Bessel
	Jid = id*r*jid                          # Ricatti-Bessel
	h1id   = MieSeries.shankel(ll, id*r)    # Spherical Hankel of first kind
	H1id = id*r*h1id                        # Ricatti-Hankel of first kind
    djid = ((ll+1).*MieSeries.sbessel(ll-1, id*r)  - ll.*MieSeries.sbessel(ll+1, id*r))  ./ (2*ll+1)  # Derivative of Spherical Bessel
	dJid = id*r*djid                                                                                  # Derivative of Ricatti-Bessel
    dh1id = ((ll+1).*MieSeries.shankel(ll-1, id*r)  - ll.*MieSeries.shankel(ll+1, id*r))  ./ (2*ll+1) # Derivative of Spherical Hankel of first kind
	dH1id = id*r*dh1id                                                                                # Derivative of Ricatti-Hankel of first kind

	vnxX = Zn* 2.0*Jid.*H1id;
	vX = Zn* 2.0*(dJid.*dH1id + ll.* (ll+1).* jid .*h1id);
	return (vnxX, vX);
end

function exciteIBC(obj, omega, pol)

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
    jk0   = MieSeries.sbessel(ll, k0*r)
    djk0 = ((ll+1).*MieSeries.sbessel(ll-1, k0*r)  - ll.*MieSeries.sbessel(ll+1, k0*r))  ./ (2*ll+1)
    h1k0   = MieSeries.shankel(ll, k0*r)
    dh1k0 = ((ll+1).*MieSeries.shankel(ll-1, k0*r)  - ll.*MieSeries.shankel(ll+1, k0*r))  ./ (2*ll+1)
	h2k0 = conj(h1k0);      # Spherical Hankel of second kind
	dh2k0 = conj(dh1k0);    # Derivative of Spherical Hankel of second kind

	(vnxX, vX) = relationSH(obj, ll);
	
    # compute the expansion coefficients for all requested polarisations
    for p in pols
        for l in ll
            a1 = -0.5*(-im)^l * sqrt(4*pi*(2*l+1)) * (-1-(h2k0[l]-im*vnxX[l]*dh2k0[l])/(h1k0[l]-im*vnxX[l]*dh1k0[l]))
            a2 = -0.5*p*(-im)^l * sqrt(4*pi*(2*l+1)) * (-1-(vX[l]*h2k0[l]-im*dh2k0[l])/(vX[l]*h1k0[l]-im*dh1k0[l]))
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

    MieSeries.Expansion(K, omega, la1, la2, ra1, ra2)
end

function MieSeries.excite(obj :: NdeltaSphere, omega, pol)
	exciteIBC(obj, omega, pol);
end

function MieSeries.excite(obj :: IBCSphere, omega, pol)
	exciteIBC(obj, omega, pol);
end

function plotRCS(sphere, omega, ls, col, mk, lb)
	pol = :x
	exc = excite(sphere, omega, pol)

	thetaArray = linspace(0.0, pi, 180)
	phiArray = linspace(0.0, 0.0, 1)
	P = [ [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)] for theta in thetaArray for phi in phiArray]

	N = length(P)
	E = zeros(N)
	H = zeros(N)

	for (i,p) in enumerate(P)
		Ex, Ey, Ez, Hx, Hy, Hz = MieSeries.farfield(sphere, exc, p)
		E[i] = norm([Ex,Ey,Ez])
		H[i] = norm([Ex,Ey,Ez])
	end

	RCS = 20*log10(sqrt(4.0*pi)*abs(E));
	plot(180.0/pi*thetaArray, RCS, linestyle=ls, color=col, marker=mk, label=lb);
end

######################################################################

const epsilon0 = 8.8541878176e-12;
const mu0 = pi*4e-7;
const eta0 = sqrt(mu0/epsilon0);
const c0 = 1.0/sqrt(epsilon0*mu0);

wavelength = 2.0;
frequency = c0/wavelength;
omega = 2.0 * pi * frequency;
k = omega * sqrt(epsilon0*mu0);

R = 1.0;
delta = 0.3;
Zs = 0.1-0.1*im; # Z divided by eta0

sphereNdelta = NdeltaSphere((Zs), delta, R, epsilon0, mu0)
sphereIBC = IBCSphere((Zs), R, epsilon0, mu0)
spherePEC = PECSphere(R, epsilon0, mu0)

using PyPlot
figure();
plotRCS(spherePEC, omega, "-", "green", "^", "RCS: PEC")
plotRCS(sphereNdelta, omega, "-", "blue", "o", "RCS: IBC Ndelta")
plotRCS(sphereIBC, omega, "-", "red", "x", "RCS: IBC nx")
legend(loc="best");

