using HallThruster: HallThruster as het
using LinearAlgebra: LinearAlgebra as linalg

A = het.Tridiagonal(ones(3), -2.6 * ones(4), ones(3))
A1 = linalg.Tridiagonal(copy(A.dl), copy(A.d), copy(A.du))
b = [-240.0, 0, 0, -150]
b1 = copy(b)
y = similar(b)
het.tridiagonal_solve!(y, A, b)
@test A1 \ b1 == y
