using DataAssim
using LinearAlgebra
using PyPlot
using Random
using Test

function check(ℳ::AbstractModel,n,t = 0,ϵ = 1e-5)
    dx = randn(n)
    x = randn(n)
    dx2 = randn(n)

    @test (ℳ(t,x + ϵ*dx) - ℳ(t,x - ϵ*dx)) / (2*ϵ)  ≈ tgl(ℳ,t,x,dx) atol=10*ϵ^2
    @test dx2 ⋅ tgl(ℳ,t,x,dx) ≈ adj(ℳ,t,x,dx2) ⋅ dx   atol=1e-7

    dX = randn(n,3)
    MdX = tgl(ℳ,t,x,dX)
    @test tgl(ℳ,t,x,dX[:,1]) ≈ MdX[:,1]
end

ℳ = ModelMatrix(2*I)

x = randn(4)
@test ℳ(0,x) ≈ 2*x
@test tgl(ℳ,0,0,x) ≈ 2*x
@test adj(ℳ,0,0,x) ≈ 2*x
check(ℳ,4)


ℳ = ModelFun((t,x,η) -> 2*x,(t,x,dx) -> 2*dx,(t,x,dx) -> 2*dx)

x = randn(4)
@test ℳ(0,x) ≈ 2*x
@test tgl(ℳ,0,0,x) ≈ 2*x
@test adj(ℳ,0,0,x) ≈ 2*x
check(ℳ,4)

# state size x
n = 2;

# number of observation per time instance
m = 1;

# observation operator
H = [1 0];
𝓗 = ModelMatrix(H)

# initial condition
xi = [1; 1];

# error covariance of the initial condition
Pi = Matrix(I,n,n)

# error covariance of the observations
R = 0.1 * Matrix(I,m,m)

ℳ = ModelMatrix([1 0.1; -0.1 1])

nmax = 100;
no = 10:5:nmax
yo = randn(m,length(no))


xai, = fourDVar(xi,Pi,ℳ,yo,R,H,nmax,no);
Q = zeros(size(Pi))
xa, = FreeRun(ℳ,xai,Q,H,nmax,no)

#𝓗
#𝓜
xa3, = KalmanFilter(xi,Pi,ℳ,Q,yo,R,H,nmax,no);
# should be ~0
@test xa[:,end] ≈ xa3[:,end]  atol=1e-5
time = 1:nmax+1

plot(time,xa3[1,:],label="KF")
plot(time,xa[1,:],label="4DVar")
#plot(time[no],yo[1,:],"*";label="observations")
legend()
PyPlot.grid("on")

xt,xfree,xa,yt,yo,diag_ = TwinExperiment(ℳ,xit,Pi,Q,R,H,nmax,no,method)
