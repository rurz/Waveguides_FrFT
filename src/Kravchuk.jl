using HypergeometricFunctions
using SpecialFunctions

"""Kravchuk.jl contains the functions to evolve a set of N + 1 amplitude points
over the discrete and finite Kravchuk oscillator realization"""

# Kravchuk symmetric polynomial definition. This particular functions is using
# the resource given in HypergeometricFunctions.jl, the method 2, who avoids
# floating-points problems, but has a flaw for the definitions in the corners
# when n = 2j and m = j (energy and position respectively).
kravd(n, m, j) = HypergeometricFunctions._₂F₁general2(-n, -m, -2j, 2)

# Normalizing terms and binomial re-definition
A(n, j) = (-1)^n / 2^j
binom(a, j) = binomial(2j, a)

# Primitive Kravchuk symmetric function. Needs a fix for the flaw described above
kravp(n, m, j) = A(n, j) * sqrt(binom(n, j) * binom(m + j, j)) * kravd(n, m + j, j)

# Conditional ad-hoc solution for the problem in the _2F1 function
function krav(n, m, j)
    if n == 2j && m == j
        return kravp(0, -j, j)
    else
        return kravp(n, m, j)
    end
end

# Exponential counting the rotation angles of the fractional evolution
kerk(n, α) = exp(-1im * pi * α * n / 2)

# Evolution kernel, this is intended to be used as standalone
frk(n, m, α, j) = sum([krav(l, n, j) * kerk(l, α) * krav(l, m, j) for l in 0:2j])

# Evolution of initial amplitudes in the vector v. This functions is intended
# to be used when the evolution is required, just given v of size 2j + 1, and
# the numbers of steps p
function evol(v, p, j)
    θ = range(0, stop = 2, length = p)
    em = zeros(Float64, (2j + 1, p))
    for r in 1:p
        for s in -j:j
        em[s + 1 + j, r] = abs2(sum([frk(s, m - j - 1, θ[r], j) * v[m] for m in 1:(2j +1)]))
        end
    end
    return em
end
