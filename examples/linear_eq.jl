# Examples from Luc Jaulin, Michel Kieffer, Olivier Didrit and Eric Walter - Applied Interval Analysis

using IntervalArithmetic, IntervalRootFinding, StaticArrays

A = [4..5 -1..1 1.5..2.5; -0.5..0.5 -7.. -5 1..2; -1.5.. -0.5 -0.7.. -0.5 2..3]
sA = SMatrix{3}{3}(A)
mA = MMatrix{3}{3}(A)

b = [3..4, 0..2, 3..4]
sb = SVector{3}(b)
mb = MVector{3}(b)

p = fill(-1e16..1e16, 3)

rts = gauss_seidel_interval!(p, A, b, precondition=true) # Gauss-Seidel Method; precondition=true by default
rts = gauss_seidel_interval!(p, sA, sb, precondition=true) # Gauss-Seidel Method; precondition=true by default
rts = gauss_seidel_interval!(p, mA, mb, precondition=true) # Gauss-Seidel Method; precondition=true by default

rts = gauss_seidel_interval(A, b, precondition=true) # Gauss-Seidel Method; precondition=true by default
rts = gauss_seidel_interval(sA, sb, precondition=true) # Gauss-Seidel Method; precondition=true by default
rts = gauss_seidel_interval(mA, mb, precondition=true) # Gauss-Seidel Method; precondition=true by default

rts = gauss_seidel_contractor!(p, A, b, precondition=true) # Gauss-Seidel Method (Vectorized); precondition=true by default
rts = gauss_seidel_contractor!(p, sA, sb, precondition=true) # Gauss-Seidel Method (Vectorized); precondition=true by default
rts = gauss_seidel_contractor!(p, mA, mb, precondition=true) # Gauss-Seidel Method (Vectorized); precondition=true by default

rts = gauss_seidel_contractor(A, b, precondition=true) # Gauss-Seidel Method (Vectorized); precondition=true by default
rts = gauss_seidel_contractor(sA, sb, precondition=true) # Gauss-Seidel Method (Vectorized); precondition=true by default
rts = gauss_seidel_contractor(mA, mb, precondition=true) # Gauss-Seidel Method (Vectorized); precondition=true by default

rts = gauss_elimination_interval!(p, A, b, precondition=true) # Gaussian Elimination; precondition=true by default
rts = gauss_elimination_interval(A, b, precondition=true) # Gaussian Elimination; precondition=true by default
