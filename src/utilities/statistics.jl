mean(x) = sum(x) / length(x)

function var(x)
    μ = mean(x)
    return mean((_x - μ)^2 for _x in x)
end

std(x) = sqrt(var(x))
