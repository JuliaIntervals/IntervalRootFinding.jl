using IntervalArithmetic
using IntervalRootFinding
using StaticArrays

# 10 dimensional problem from https://discourse.julialang.org/t/solvers-fail-on-nonlinear-problem-which-has-solutions/20051
function f10d(X)
    p = [3600, 18, 18, 3600, 3600, 3600, 1800, 3600, 18, 18, 18, 1800, 18, 36, 11, 180, 0.7, 0.4, 30, 0.2, 4, 4.5, 0.4]
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = X

    return SVector(
        -p[17] * x1 + -2 * p[1] * (x1 ^ 2 / 2) + 2 * p[2] * x2 + p[21] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x2 + p[1] * (x1 ^ 2 / 2) + -1 * p[2] * x2 + -1 * p[4] * x2 * x4 + p[9] * x3 + p[14] * x3 + -1 * p[6] * x2 * x7 + p[11] * x8,
        -p[17] * x3 + p[4] * x2 * x4 + -1 * p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[14] * x3 + p[15] * x5 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7,
        -p[17] * x4 + -1 * p[4] * x2 * x4 + p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7 + p[16] * x9 + p[22] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x5 + p[5] * x3 * x4 + -1 * p[10] * x5 + -1 * p[15] * x5,
        -p[17] * x6 + p[14] * x3 + p[15] * x5 + -1 * p[8] * x6 * x10 + p[13] * x9,
        -p[17] * x7 + -1 * p[6] * x2 * x7 + p[11] * x8 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7 + p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x8 + p[6] * x2 * x7 + -1 * p[11] * x8 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7,
        -p[17] * x9 + p[8] * x6 * x10 + -1 * p[13] * x9 + -1 * p[16] * x9,
        x10 + x9 - p[23]
    )
end

# Approximate solution
sol10d = SVector(0.11584484298713685, 
  0.9455230206055336,   
  1.4935460680512724,   
  0.014733756971650408, 
  2.6673387628301497,   
 15.827943289607886,    
  0.03777227940579977,  
  5.088785613357363,    
  0.3986099863617728,   
  0.0013900136199035832
)

# This input is close enough to the solution to converge in less than 1e6 iterations
X10d_close = SVector{10}(interval(x - 0.0002, x + 0.0002) for x in sol10d)

# This input is too large, and requires more than 1e6 iterations to converge
X10d_large = SVector{10}(fill(interval(0, 100), 10))