using MieSeries;
using BoundaryElements;
using CompScienceMeshes;
using LinearForms;
using FixedSizeArrays;

	#################################################
	#  Xlm and Ulm testing

struct VectorSphericalHarmonicX <: BoundaryElements.Functional
	l::Int
	m::Int
end

struct VectorSphericalHarmonicU <: BoundaryElements.Functional
	l::Int
	m::Int
end

function (vshX::VectorSphericalHarmonicX)(p::MeshPointNM)
	r = cartesian(p);
	l = vshX.l;
	m = vshX.m;
	lm = l*(l+1)+m+1;
	Xx, Xy, Xz, Yx, Yy, Yz, Ux, Uy, Uz = MieSeries.Xlm(l+1, r);
	return Vec(Xx[lm], Xy[lm], Xz[lm]);
end

function (vshX::VectorSphericalHarmonicU)(p::MeshPointNM)
	r = cartesian(p);
	l = vshX.l;
	m = vshX.m;
	lm = l*(l+1)+m+1;
	Xx, Xy, Xz, Yx, Yy, Yz, Ux, Uy, Uz = MieSeries.Xlm(l+1, r);
	return Vec(Ux[lm], Uy[lm], Uz[lm]);
end

function BoundaryElements.integrand(vshX::VectorSphericalHarmonicX, tv, fv)
	return dot(tv[1], fv);
end

function BoundaryElements.integrand(vshU::VectorSphericalHarmonicU, tv, fv)
	return dot(tv[1], fv);
end

	#########################################
	# Simulation parameters

const epsilon0 = 8.8541878176e-12;
const mu0 = pi*4e-7;
const eta0 = sqrt(mu0/epsilon0);
const c0 = 1.0/sqrt(epsilon0*mu0);

wavelength = 5.0;
frequency = c0/wavelength;
omega = 2.0 * pi * frequency;
k = omega * sqrt(epsilon0*mu0);

Z = (0.7+0.3*im)*eta0;
delta = 1.2;

#l = 2;
l = 1
m = -1;

	####################################
	# Operators

R = 4.0;
h = 1.2;
gamma = meshsphere(R, h);
f = raviartthomas(gamma);
g = buffachristiansen(gamma);

#Ts = MWWeaklySingular(k);
#Th = MWHyperSingular(k);
Ts = MWWeaklySingular(im*k);
Th = MWHyperSingular(im*k);
K = MWDoubleLayer3D(im*k);
N = NCross();
I = Identity();
Tdelta = MWWeaklySingular(1/delta);

println("Ts_nxff"); Ts_nxff = 1.0*assemble(Ts,f,f);
println("Th_nxff"); Th_nxff = -1.0*assemble(Th,f,f); # T = -ikTs + Th/ik is preferred instead of -ikTs - Th/ik
println("K_nxff"); K_nxff = 1.0*assemble(K,f,f);
println("K_nxfg"); K_nxfg = assemble(K,f,g);
println("Gmix = I_nxfg"); Gmix = -assemble(N,f,g); # <nxf,g> = -<f,Ng>
println("IGmix"); IGmix = inv(Gmix);
println("I_ff = Gram"); I_ff = assemble(I,f,f);
println("IGram"); IGram = inv(I_ff);
println("N_ff = I_fnxf"); N_ff = assemble(N,f,f);
println("Nd_nxff"); Nd_nxff = (2.0/delta)*assemble(Tdelta,f,f);
# =#

	###########################################
	# Vector Spherical harmonics testing

VSH_Xlm = VectorSphericalHarmonicX(l, m);
VSH_Ulm = VectorSphericalHarmonicU(l, m);

println("Xlm_f"); Xlm_f = assemble(VSH_Xlm, f);
println("Ulm_f"); Ulm_f = assemble(VSH_Ulm, f);
println("Xlm_nxf"); Xlm_nxf = -Ulm_f;  # <nxf,X> = -<f,U>
println("Ulm_nxf"); Ulm_nxf = Xlm_f;   # <nxf,U> = <f,X>

	##########################################
	# Special functions

#   input :
# order  integer order of the functions
# arg    complex argument of the functions
#   output :
# j      Spherical Bessel
# h1     Spherical Hankel of 1st kind
# h2     Spherical Hankel of 2nd kind
# dj     Derivative of Spherical Bessel
# dh1    Derivative of Spherical Hankel of 1st kind
# dh2    Derivative of Spherical Hankel of 2nd kind
# J      Ricatti Bessel
# H1     Ricatti Hankel of 1st kind
# H2     Ricatti Hankel of 2nd kind
# dJ     Derivative of Ricatti Bessel
# dH1    Derivative of Ricatti Hankel of 1st kind
# dH2    Derivative of Ricatti Hankel of 2nd kind
function BesselFamilyFunctions(order, arg)
	j = MieSeries.sbessel(order, arg);
	h1 = MieSeries.shankel(order, arg); # j+im*y
	h2 = 2.0*j-h1;  # j-im*y equals conj(h1) only when arg is real
	dj = (order*MieSeries.sbessel(order-1, arg) - (order+1)*MieSeries.sbessel(order+1, arg)) / (2*order+1);
	dh1 = (order*MieSeries.shankel(order-1, arg) - (order+1)*MieSeries.shankel(order+1, arg)) / (2*order+1);
	dh2 = 2.0*dj-dh1;

	J = arg*j;
	H1 = arg*h1;
	H2 = 2.0*J-H1;
	dJ = j + arg*dj; # d(z*j(z))/dz
	dH1 = h1 + arg*dh1;
	dH2 = 2.0*dJ-dH1;

	return (j, h1, h2, dj, dh1, dh2, J, H1, H2, dJ, dH1, dH2);
end

(j, h1, h2, dj, dh1, dh2,
 J, H1, H2, dJ, dH1, dH2) = BesselFamilyFunctions(l, k*R);

Rid = R/(im*delta);
(j_Rid, h1_Rid, h2_Rid, dj_Rid, dh1_Rid, dh2_Rid,
 J_Rid, H1_Rid, H2_Rid, dJ_Rid, dH1_Rid, dH2_Rid) = BesselFamilyFunctions(l, Rid);

println("im*J*dH2-im*dJ*H2 = ($(im*J*dH2-im*dJ*H2)) should be 1");
println("im*J_Rid*dH2_Rid-im*dJ_Rid*H2_Rid = ($(im*J_Rid*dH2_Rid-im*dJ_Rid*H2_Rid)) should be 1");

	#########################################
	# Numerical/Analytical comparisons

# N

num_NX = N_ff*IGram*Xlm_f;
anal_NX = Ulm_f;
re_NX = norm(num_NX-anal_NX)./norm(anal_NX);
println("re_NX : $re_NX");

num_NU = N_ff*IGram*Ulm_f;
anal_NU = -Xlm_f;
re_NU = norm(num_NU-anal_NU)./norm(anal_NU);
println("re_NU : $re_NU");

# T

T_nxff = -im*k*Ts_nxff + 1.0/(im*k)*Th_nxff;

num_TX = T_nxff*IGram*Xlm_f;
anal_TX = -J*H2*Ulm_nxf;
re_TX = norm(num_TX-anal_TX)./norm(anal_TX);
println("re_TX : $re_TX");

num_TU = T_nxff*IGram*Ulm_f;
anal_TU = dJ*dH2*Xlm_nxf;
re_TU = norm(num_TU-anal_TU)./norm(anal_TU);
println("re_TU : $re_TU");

# Ts

num_TsX = Ts_nxff*IGram*Xlm_f;
anal_TsX = 1.0/(im*k)*J*H2*Ulm_nxf;
re_TsX = norm(num_TsX-anal_TsX)./norm(anal_TsX);
println("re_TsX : $re_TsX");

num_TsU = Ts_nxff*IGram*Ulm_f;
anal_TsU = 1.0/(-im*k)*(dJ*dH2+l*(l+1)*j*h2)*Xlm_nxf;
re_TsU = norm(num_TsU-anal_TsU)./norm(anal_TsU);
println("re_TsU : $re_TsU");

# Th

# Print the absolute error for the numerical Th(X) because it equals 0 analytically
# Compare it relatively to T(X)
num_ThX = Th_nxff*IGram*Xlm_f;
anal_ThX = 0*Ulm_nxf;
abserror_ThX = norm(num_ThX);
reToTX_ThX = abserror_ThX / norm(anal_TX);
println("abserror_ThX = $(abserror_ThX) (norm(anal_TX) = $(norm(anal_TX)), reToTX_ThX = $reToTX_ThX)");

num_ThU = Th_nxff*IGram*Ulm_f;
anal_ThU = -im*k*l*(l+1)*j*h2*Xlm_nxf;
re_ThU = norm(num_ThU-anal_ThU)./norm(anal_ThU);
println("re_ThU : $re_ThU");

# K

num_KX = K_nxff*IGram*Xlm_f;
anal_KX = (im*J*dH2-0.5)*Xlm_nxf;
re_KX = norm(num_KX-anal_KX)./norm(anal_KX);
println("re_KX : $re_KX");

num_KU = K_nxff*IGram*Ulm_f;
anal_KU = (-im*dJ*H2-0.5)*Ulm_nxf;
re_KU = norm(num_KU-anal_KU)./norm(anal_KU);
println("re_KU : $re_KU");

# (I/2-K)N

Nover2minusKN_nxff = (0.5*Gmix-K_nxfg)*IGmix*I_ff;

num_Nover2minusKNX = Nover2minusKN_nxff*IGram*Xlm_f;
anal_Nover2minusKNX = im*J*dH2*Ulm_nxf;
re_Nover2minusKNX = norm(num_Nover2minusKNX-anal_Nover2minusKNX)./norm(anal_Nover2minusKNX);
println("re_Nover2minusKNX : $re_Nover2minusKNX");

num_Nover2minusKNU = Nover2minusKN_nxff*IGram*Ulm_f;
anal_Nover2minusKNU = im*dJ*H2*Xlm_nxf;
re_Nover2minusKNU = norm(num_Nover2minusKNU-anal_Nover2minusKNU)./norm(anal_Nover2minusKNU);
println("re_Nover2minusKNU : $re_Nover2minusKNU");

# S = eta * T - Z * (I/2-K)N

S_nxff = eta0 * T_nxff - Z * Nover2minusKN_nxff;

num_SX = S_nxff*IGram*Xlm_f;
anal_SX = (-eta0*H2-im*Z*dH2)*J*Ulm_nxf;
re_SX = norm(num_SX-anal_SX)./norm(anal_SX);
println("re_SX : $re_SX");

num_SU = S_nxff*IGram*Ulm_f;
anal_SU = (eta0*dH2-im*Z*H2)*dJ*Xlm_nxf;
re_SU = norm(num_SU-anal_SU)./norm(anal_SU);
println("re_SU : $re_SU");

# N delta

num_NdX = Nd_nxff*IGram*Xlm_f;
anal_NdX = 2.0*J_Rid*H2_Rid*Ulm_nxf;
re_NdX = norm(num_NdX-anal_NdX)./norm(anal_NdX);
println("re_NdX : $re_NdX");

num_NdU = Nd_nxff*IGram*Ulm_f;
anal_NdU = -2.0*(dJ_Rid*dH2_Rid+l*(l+1)*j_Rid*h2_Rid)*Xlm_nxf;
re_NdU = norm(num_NdU-anal_NdU)./norm(anal_NdU);
println("re_NdU : $re_NdU");
