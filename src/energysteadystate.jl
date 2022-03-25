#eliminating the time derivative in energy equation
#see the resulting changes in the overall solution

#put in update values callback just like the potential solver, can then comment in or out

function solve_energyss!(U, params)
    #directly discretising the equation, conserves properties such as negative semidefinite etc...
    #add functionality for nonuniform cell size
    z_cell = params.z_cell
    fluids = params.fluids
    index = params.index
    N = length(z_cell)

    #need to make larger A and b, same as cells
    A = params.cache.A₁
    b = params.cache.b₁

    #Δz = z_cell[i+1] - z_cell[i]
    Δz = z_cell[3] - z_cell[2]

    mi = params.propellant.m
    νϵ = 1e7 * smooth_if(params.z_cell[end-1], params.L_ch, params.νϵ[1], params.νϵ[2], 10)
    UU = 20.0
    W = νϵ * exp(-UU / U[index.Tev, end-1]) #without Tev multiplication here
    νϵ_tot = #=U[1, end]/mi*params.loss_coeff(U[index.Tev, end])=# + W
    ue = U[index.ue, end-1]

    A.d[1] = 1.0
    A.du[1] = 0.0
    A.d[N] = 5/3 /Δz - (νϵ_tot) / ue
    A.dl[N-1] = -5/3 /Δz

    b[1] = params.Te_L
    b[end] = -U[index.grad_ϕ, end-1]

    mi = m(fluids[1])

    #=@turbo=# for i in 2:(N-1)

        ue = U[index.ue, i]
        
        if abs(ue) < 100 
            if ue < 0 
                ue = -100.0
            else 
                ue = 100.0
            end
        end

        #get loss term for steady state energy equation
        mi = params.propellant.m
        νϵ = 1e7 * smooth_if(params.z_cell[i], params.L_ch, params.νϵ[1], params.νϵ[2], 10)
        UU = 20.0
        W = νϵ * exp(-UU / U[index.Tev, i]) #without Tev multiplication here
        νϵ_tot = #=U[1, i]/mi*params.loss_coeff(U[index.Tev, i])=# + W
        

        #direct discretization, h to each side
        A.dl[i - 1] = -5/6 / Δz
        A.d[i] = -(νϵ_tot) / ue
        A.du[i] = 5/6 / Δz

        #=
        @show(W)
        @show(U[1, i]/mi*params.loss_coeff(U[index.Tev, i]))
        @show(A.d[i])=#

        #source term, h to each side
        b[i] = -U[index.grad_ϕ, i]
    end

    U[index.Tev, :] = inv(A)*b
    #@show(U[index.Tev, :])

    for i in 1:N
        U[index.Tev, i] = max(3.0, U[index.Tev, i])
    end

    return U[index.Tev, :] .* U[index.ne, :]
end