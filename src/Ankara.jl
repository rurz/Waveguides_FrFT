using SpecialMatrices
using LinearAlgebra
using SparseArrays

"""Ankara.jl contain the functions to develop the description of the
discrete and finite harmonic oscillator in terms of the Harper equations."""

# Δ1(j) gives the circulant matrix of size d_{j} × d_{j} where d_{j} = 2j + 1
function Δ1(j)
	Δ1 = zeros(Float64, 2j + 1)
	Δ1[1] = -2
	Δ1[2] = 1
	Δ1[end] = 1
	Δ1 = SpecialMatrices.Circulant(Δ1)
end

# Δ2(j) gives the diagonal matrix of size d_{j} × d_{j} whose diagonal entries
# are oscillating terms
function Δ2(j)
	Δ2 = zeros(Float64, 2j + 1)
	for i in -j:j
		Δ2[i + 1 + j] = -4 * (sin(pi * i / (2j)))^2
	end
	return Δ2 = LinearAlgebra.diagm(Δ2)
end

# If necessary, compose the Hamiltonian (Harper's equation) H from who the
# eigenvectors and eigenvalues are obtained
H(j) = -(1/2.0)*(Δ1(j) + Δ2(j))

# Eigenvectors of H
h(j) = LinearAlgebra.eigvecs(H(j));

# Eigenvalues of H. Altough, this eigvenvalues are not strictly used because
# of the symmetry importation below
k(j) = LinearAlgebra.eigvals(H(j));

# Lineal monotonic increasing values who serves as the (new) eigenvalues
# in the symmetry importation to construct the evolution operator
η(j) = [l for l in 0:(2j + 1)];

# Evolution operator
G(j, m, μ, t) = sum([h(j)[m + j + 1, l + 1] * exp(-1im * pi * t * η(j)[l + 1]/2.0) * h(j)[μ + j + 1, l + 1] for l in 0:2j]);

# If necessary, this function gives the matrix of the evolution operator
# in function of t. IT IS REALLY SLOW
P(j, t) = reshape([G(j, m, μ, t) for m in -j:j for μ in -j:j], (2j + 1, 2j + 1))
