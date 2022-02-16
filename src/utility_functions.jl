"""
    smooth_max(x, y, k=10)
Computes a smooth approximation to max(x, y)
"""
function smooth_max(x,y,k = 10)
    (x*exp(k*x) + y*exp(k*y)) / (exp(k*x) + exp(k*y))
end

"""
    smooth_min(x, y, k=10)
Compute a smooth approximation to min(x, y)
"""
smooth_min(x,y,k=10) = smooth_max(x, y, -k)

"""
    smooth_if_gt(x, cutoff, v1, v2, k=10)
Computes an analytic approximation to x < cutoff ? v1 : v2
"""
smooth_if(x, cutoff, v1, v2, k = 10) = 0.5*((v2-v1)*tanh(k*(x-cutoff)) + v1+v2)
