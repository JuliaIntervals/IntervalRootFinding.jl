# examples from
# M. W. Smiley and C. Chun, J. Comput. Appl. Math. **137**, 293 (2001).
# https://doi.org/10.1016/S0377-0427(00)00711-1

module SmileyExample22

using IntervalArithmetic
using StaticArrays

const title = "Smiley and Chun (2001), Example 2.2"

f(x) = SVector(
    x[1]^2 + 4 * x[2]^2 - 4,
    x[2] * (x[1] - 1.995) * (x[2] - x[1]^2) * (x[2] - x[1] + 1)
)

# contains all 8 reported roots
const region = IntervalBox(-3..3, -3..3)

# only three roots are reported explicitely:
# (2, 0) and (1.995, ±0.071)

end

module SmileyExample52

using IntervalArithmetic
using StaticArrays

const title = "Smiley and Chun (2001), Example 5.2"

const c0  =  1
const a0  =  0.5
const b0  =  0.5
const d01 =  0.2
const d02 = -0.7
const r1  =  1
const r2  =  2
const m   =  3
const θs  = [i * pi / m for i in 1:m]
_v(θ) = -(a0 * cos(θ) + b0 * sin(θ)) / c0
const vs = _v.(θs)

_a(θ, v) = b0 * v - c0 * sin(θ)
const as = _a.(θs, vs)
_b(θ, v) = c0 * cos(θ) - a0 * v
const bs = _b.(θs, vs)
_c(θ) = a0 * sin(θ) - b0 * cos(θ)
const cs = _c.(θs)
_d(c) = d01 * c / c0
const ds = _d.(cs)

@show as, bs, cs, ds

f(x) = SVector(
    (x[1]^2 + x[2]^2 + x[3]^2 - r1^2) * (x[1]^2 + x[2]^2 + x[3]^2 - r2^2),
    (a0 * x[1] + b0 * x[2] + c0 * x[3] - d01) *
        (a0 * x[1] + b0 * x[2] + c0 * x[3] - d02),
    prod(as[i] * x[1] + bs[i] * x[2] + cs[i] * x[3] - ds[i] for i in 1:m)
)

# contains all 24 reported roots
const region = IntervalBox(-3..3, -3..3, -3..3)

const known_roots = [
    IntervalBox(-1.933009 ± 1e-6, -0.300000 ± 1e-6,  0.416504 ± 1e-6),
    IntervalBox(-1.701684 ± 1e-6,  0.000000 ± 1e-6,  1.050842 ± 1e-6),
    IntervalBox(-1.258803 ± 1e-6,  1.360696 ± 1e-6, -0.750946 ± 1e-6),
    IntervalBox(-1.044691 ± 1e-6, -1.589843 ± 1e-6,  0.617267 ± 1e-6),
    IntervalBox(-0.996600 ± 1e-6,  1.726162 ± 1e-6, -0.164780 ± 1e-6),

    IntervalBox(-0.951026 ± 1e-6, -0.300000 ± 1e-6, -0.074486 ± 1e-6),
    IntervalBox(-0.800000 ± 1e-6, -0.000000 ± 1e-6,  0.600000 ± 1e-6),
    IntervalBox(-0.776373 ± 1e-6, -1.344718 ± 1e-6,  1.260546 ± 1e-6),
    IntervalBox(-0.717665 ± 1e-6,  0.423418 ± 1e-6, -0.552876 ± 1e-6),
    IntervalBox(-0.592072 ± 1e-6, -0.805884 ± 1e-6, -0.001021 ± 1e-6),

    IntervalBox(-0.499927 ± 1e-6,  0.865900 ± 1e-6,  0.017013 ± 1e-6),
    IntervalBox(-0.360640 ± 1e-6, -0.624646 ± 1e-6,  0.692643 ± 1e-6),
    IntervalBox( 0.082249 ± 1e-6, -0.962075 ± 1e-6, -0.260086 ± 1e-6),
    IntervalBox( 0.085220 ± 1e-6,  0.367221 ± 1e-6, -0.926221 ± 1e-6),
    IntervalBox( 0.453788 ± 1e-6,  0.785984 ± 1e-6, -0.419886 ± 1e-6),

    IntervalBox( 0.464511 ± 1e-6, -0.804557 ± 1e-6,  0.370022 ± 1e-6),
    IntervalBox( 0.511026 ± 1e-6, -0.300000 ± 1e-6, -0.805513 ± 1e-6),
    # the following two roots are suspect, first column probably reported in error
    #IntervalBox( 0.623386 ± 1e-6,  1.151180 ± 1e-6, -1.544510 ± 1e-6),
    #IntervalBox( 0.869521 ± 1e-6, -1.899353 ± 1e-6, -0.062016 ± 1e-6),
    IntervalBox( 0.537839 ± 1e-6,  1.151180 ± 1e-6, -1.544510 ± 1e-6),
    IntervalBox( 0.623386 ± 1e-6, -1.899353 ± 1e-6, -0.062016 ± 1e-6),
    IntervalBox( 0.869521 ± 1e-6,  1.506056 ± 1e-6, -0.987788 ± 1e-6),

    IntervalBox( 0.960000 ± 1e-6,  0.000000 ± 1e-6, -0.280000 ± 1e-6),
    IntervalBox( 0.961183 ± 1e-6, -1.664819 ± 1e-6,  0.551817 ± 1e-6),
    IntervalBox( 1.493009 ± 1e-6, -0.300000 ± 1e-6, -1.296504 ± 1e-6),
    IntervalBox( 1.861684 ± 1e-6,  0.000000 ± 1e-6, -0.730842 ± 1e-6),
]

end

# example 5.4, rescaled form
module SmileyExample54

using IntervalArithmetic
using StaticArrays

const title = "Smiley and Chun (2001), Example 5.4"

const c11 =  1.069e-5
const c12 =  2e2
const c13 =  1e5
const c14 = -1.8e5
const c15 = -1.283e-4
const c21 =  2e-2
const c22 =  1e1
const c23 = -1e1

f(t) = SVector(
    c11 * t[1]^4 + c12 * t[1]^3 * t[2] + c13 * t[1]^3 + c14 * t[1] + c15,
    c21 * t[1] * t[2]^2 + c22 * t[2]^2 + c23
    )

# contains all 7 reported roots
const region = IntervalBox(-5.1e2..1.4, -5.1e2..1.1)

const known_roots = [
    IntervalBox(-1.34298 ± 1e-5, -1.00134 ± 1e-5),
    IntervalBox(-1.34030 ± 1e-5,  1.00134 ± 1e-5),
    IntervalBox( 1.34030 ± 1e-5,  0.99866 ± 1e-5),
    IntervalBox( 1.34298 ± 1e-5, -0.99866 ± 1e-5),
    # following three roots are not reported precisely
    IntervalBox(-7e-10 ± 1e-10,  1 ± 1e-5),
    IntervalBox(-7e-10 ± 1e-10, -1 ± 1e-5),
    IntervalBox(-5e2 ± 1, -5e2 ± 1)
]

end

module SmileyExample55

using IntervalArithmetic
using StaticArrays

const title = "Smiley and Chun (2001), Example 5.5"

const μ1 = pi / 10
const μ2 = pi / 5
const α  = 5
const D1 = exp(-2μ1)
const D2 = exp(-2μ2)
const C1 = (1 - D1) / 2μ1
const C2 = (1 - D2) / 2μ2

g(x) = SVector(C1 * (x[3] - α * sin(x[1]) * cos(x[2])) + x[1],
               C2 * (x[4] - α * cos(x[1]) * sin(x[2])) + x[2],
               D1 * (x[3] - α * sin(x[1]) * cos(x[2])),
               D2 * (x[4] - α * cos(x[1]) * sin(x[2])))

f(x) = (g ∘ g)(x) .- x

# contains all 41 reported roots of f
const region =
    IntervalBox(-1.02pi..1.02pi, -1.02pi..1.02pi, -0.5pi..0.5pi, -0.5pi..0.5pi)

# roots from
# C. S. Hsu and R. S. Guttalu, Trans. ASME **50**, 858 (1983).
# http://dx.doi.org/10.1115/1.3167157
const known_roots = [
    # period 1 points
    IntervalBox(-pi    ± 0, -pi    ± 0, 0 ± 0, 0 ± 0),
    IntervalBox(-pi    ± 0,  0     ± 0, 0 ± 0, 0 ± 0),
    IntervalBox(-pi    ± 0,  pi    ± 0, 0 ± 0, 0 ± 0),
    IntervalBox(-0.5pi ± 0, -0.5pi ± 0, 0 ± 0, 0 ± 0),
    IntervalBox(-0.5pi ± 0,  0.5pi ± 0, 0 ± 0, 0 ± 0),

    IntervalBox( 0     ± 0, -pi    ± 0, 0 ± 0, 0 ± 0),
    IntervalBox( 0     ± 0,  0     ± 0, 0 ± 0, 0 ± 0),
    IntervalBox( 0     ± 0,  pi    ± 0, 0 ± 0, 0 ± 0),
    IntervalBox( 0.5pi ± 0, -0.5pi ± 0, 0 ± 0, 0 ± 0),
    IntervalBox( 0.5pi ± 0,  0.5pi ± 0, 0 ± 0, 0 ± 0),

    IntervalBox( pi    ± 0, -pi    ± 0, 0 ± 0, 0 ± 0),
    IntervalBox( pi    ± 0,  0     ± 0, 0 ± 0, 0 ± 0),
    IntervalBox( pi    ± 0,  pi    ± 0, 0 ± 0, 0 ± 0),
    # period 2 points
    IntervalBox( (0.33419 ± 1e-5)pi,  0 ± 0,  (0.48025 ± 1e-5)pi, 0 ± 0),
    IntervalBox(-(0.33419 ± 1e-5)pi,  0 ± 0, -(0.48025 ± 1e-5)pi, 0 ± 0),

    IntervalBox( 0 ± 0,  (0.24702 ± 1e-5)pi, 0 ± 0,  (0.24699 ± 1e-5)pi),
    IntervalBox( 0 ± 0, -(0.24702 ± 1e-5)pi, 0 ± 0, -(0.24699 ± 1e-5)pi),
    IntervalBox( (0.09648 ± 1e-5)pi,  (0.18318 ± 1e-5)pi,  (0.13864 ± 1e-5)pi,  (0.18316 ± 1e-5)pi),
    IntervalBox(-(0.09648 ± 1e-5)pi, -(0.18318 ± 1e-5)pi, -(0.13864 ± 1e-5)pi, -(0.18316 ± 1e-5)pi),
    IntervalBox(-(0.09648 ± 1e-5)pi,  (0.18318 ± 1e-5)pi, -(0.13864 ± 1e-5)pi,  (0.18316 ± 1e-5)pi),

    IntervalBox( (0.09648 ± 1e-5)pi, -(0.18318 ± 1e-5)pi,  (0.13864 ± 1e-5)pi, -(0.18316 ± 1e-5)pi),
    IntervalBox(-(0.66580 ± 1e-5)pi, -(1 ± 0)pi,  (0.48025 ± 1e-5)pi, 0 ± 0),
    IntervalBox( (0.66580 ± 1e-5)pi, -(1 ± 0)pi, -(0.48025 ± 1e-5)pi, 0 ± 0),
    IntervalBox(-(0.66580 ± 1e-5)pi,  (1 ± 0)pi,  (0.48025 ± 1e-5)pi, 0 ± 0),
    IntervalBox( (0.66580 ± 1e-5)pi,  (1 ± 0)pi, -(0.48025 ± 1e-5)pi, 0 ± 0),

    IntervalBox( (1 ± 0)pi,  -(0.75298 ± 1e-5)pi, 0 ± 0,  (0.24699 ± 1e-5)pi),
    IntervalBox( (1 ± 0)pi,   (0.75298 ± 1e-5)pi, 0 ± 0, -(0.24699 ± 1e-5)pi),
    IntervalBox(-(1 ± 0)pi,  -(0.75298 ± 1e-5)pi, 0 ± 0,  (0.24699 ± 1e-5)pi),
    IntervalBox(-(1 ± 0)pi,   (0.75298 ± 1e-5)pi, 0 ± 0, -(0.24699 ± 1e-5)pi),
    IntervalBox( (0.90352 ± 1e-5)pi,  (0.81682 ± 1e-5)pi, -(0.13864 ± 1e-5)pi, -(0.18316 ± 1e-5)pi),

    IntervalBox(-(0.90352 ± 1e-5)pi, -(0.81682 ± 1e-5)pi,  (0.13864 ± 1e-5)pi,  (0.18316 ± 1e-5)pi),
    IntervalBox( (0.90352 ± 1e-5)pi, -(0.81682 ± 1e-5)pi, -(0.13864 ± 1e-5)pi,  (0.18316 ± 1e-5)pi),
    IntervalBox(-(0.90352 ± 1e-5)pi,  (0.81682 ± 1e-5)pi,  (0.13864 ± 1e-5)pi, -(0.18316 ± 1e-5)pi),
    IntervalBox( (0.34989 ± 1e-5)pi, -(0.64408 ± 1e-5)pi, -(0.21572 ± 1e-5)pi, -(0.14406 ± 1e-5)pi),
    IntervalBox( (0.65011 ± 1e-5)pi, -(0.35592 ± 1e-5)pi,  (0.21572 ± 1e-5)pi,  (0.14406 ± 1e-5)pi),

    IntervalBox( (0.34989 ± 1e-5)pi,  (0.64408 ± 1e-5)pi, -(0.21572 ± 1e-5)pi,  (0.14406 ± 1e-5)pi),
    IntervalBox( (0.65011 ± 1e-5)pi,  (0.35592 ± 1e-5)pi,  (0.21572 ± 1e-5)pi, -(0.14406 ± 1e-5)pi),
    IntervalBox(-(0.34989 ± 1e-5)pi, -(0.64408 ± 1e-5)pi,  (0.21572 ± 1e-5)pi, -(0.14406 ± 1e-5)pi),
    IntervalBox(-(0.65011 ± 1e-5)pi, -(0.35592 ± 1e-5)pi, -(0.21572 ± 1e-5)pi,  (0.14406 ± 1e-5)pi),
    IntervalBox(-(0.34989 ± 1e-5)pi,  (0.64408 ± 1e-5)pi,  (0.21572 ± 1e-5)pi,  (0.14406 ± 1e-5)pi),

    IntervalBox(-(0.65011 ± 1e-5)pi,  (0.35592 ± 1e-5)pi, -(0.21572 ± 1e-5)pi, -(0.14406 ± 1e-5)pi)
]

end
