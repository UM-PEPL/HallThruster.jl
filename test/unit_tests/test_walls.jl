@testset "Wall loss tests" begin
    no_losses = HallThruster.NoWallLosses()

    mi = HallThruster.Xenon.m
    me = HallThruster.me

    αin = 0.3
    αout = 1.0

    landmark_losses = HallThruster.ConstantSheathPotential(20.0, αin, αout)

    Tev = 4.0
    ne = 1e18
    cache = (;
        ne = [ne, ne, ne, ne],
        Tev = [Tev, Tev, Tev, Tev],
        Z_eff = [1.0, 1.0, 1.0, 1.0],
        ni = [ne ne ne ne],
        γ_SEE = [0.0, 0.0, 0.0, 0.0],
        νew = [0.0, 0.0, 0.0, 0.0],
    )
    transition_function = HallThruster.StepFunction()
    config = (;
        thruster = HallThruster.SPT_100,
        transition_function,
        propellant = HallThruster.Xenon, ncharge = 1,
        electron_plume_loss_scale = 0.0
    )
    index = (;ρi = [1], nϵ = 2)

    L_ch = HallThruster.geometry_SPT_100.channel_length
    A_ch = HallThruster.channel_area(HallThruster.SPT_100)

    grid = HallThruster.generate_grid(HallThruster.geometry_SPT_100, 2, (0, 2 * L_ch))

    z_cell = grid.cell_centers
    z_edge = grid.edges

    # Fill up cell lengths and magnetic field vectors
    Δz_cell = zeros(length(z_cell))
    Δz_edge = zeros(length(z_edge))
    for (i, z) in enumerate(grid.cell_centers)
        if firstindex(z_cell) < i < lastindex(z_cell)
            Δz_cell[i] = z_edge[HallThruster.right_edge(i)] - z_edge[HallThruster.left_edge(i)]
        elseif i == firstindex(z_cell)
            Δz_cell[i] = z_edge[begin+1] - z_edge[begin]
        elseif i == lastindex(z_cell)
            Δz_cell[i] = z_edge[end] - z_edge[end-1]
        end
    end

    for i in eachindex(z_edge)
        Δz_edge[i] = z_cell[i+1] - z_cell[i]
    end


    mi_kr = HallThruster.Krypton.m
    γmax = 1 - 8.3 * sqrt(me/mi)
    γmax_kr = 1 - 8.3 * sqrt(me/mi_kr)

    params = (;
        cache, config, index, z_cell = grid.cell_centers,
        z_edge = grid.edges, L_ch = L_ch, A_ch = A_ch,
        Δz_edge, Δz_cell, γ_SEE_max = γmax
    )

    ρi = ne * mi
    nϵ = 3/2 * ne * Tev

    U = [
        ρi ρi ρi ρi
        nϵ nϵ nϵ nϵ
    ]

    @test HallThruster.wall_power_loss(no_losses, U, params, 2) == 0.0
    @test HallThruster.wall_power_loss(landmark_losses, U, params, 2) ≈ αin *  6.0 * 1e7 * exp(-20 / 6.0)
    @test HallThruster.wall_power_loss(landmark_losses, U, params, 3) ≈ αout *  6.0 * 1e7 * exp(-20 / 6.0)

    γ1 = 0.5
    γ2 = 1.0
    γsat = 100.0
    @test HallThruster.sheath_potential(Tev, γ1, mi) == Tev*log(0.5*(1-γ1)*sqrt(2*mi/π/me))
    @test HallThruster.sheath_potential(Tev, γ2, mi) == Tev*log(0.5*(1-γ2)*sqrt(2*mi/π/me))
    #-1.02 * Tev

    BN = HallThruster.BoronNitride
    @test HallThruster.SEE_yield(BN, 10.0, γmax) ≈ BN.σ₀ + 1.5 * 10.0/BN.ϵ_star * (1 - BN.σ₀)
    # Check that SEE properly obeys space charge limit at high electron temps
    @test HallThruster.SEE_yield(BN, 9000.0, γmax) ≈ γmax


    @test HallThruster.SEE_yield(BN, 9000.0, γmax_kr) ≈ γmax_kr

    γ = HallThruster.SEE_yield(BN, Tev, γmax)

    Vs = HallThruster.sheath_potential(Tev, γ, mi)

    α = 1/4
    sheath_model = HallThruster.WallSheath(BN, α)
    Δr = HallThruster.geometry_SPT_100.outer_radius - HallThruster.geometry_SPT_100.inner_radius
    Δz = params.z_edge[2] - params.z_edge[1]
    V_cell = A_ch * Δz

    @test HallThruster.freq_electron_wall(no_losses, U, params, 2) == 0.0
    @test HallThruster.freq_electron_wall(no_losses, U, params, 3) == 0.0

    @test HallThruster.freq_electron_wall(landmark_losses, U, params,  2) == 1e7
    @test HallThruster.freq_electron_wall(landmark_losses, U, params, 3) == 0.0e7

    γ = HallThruster.SEE_yield(BN, Tev, γmax)
    νiw = α * sqrt(HallThruster.e * Tev / mi) / Δr * γ
    νew = νiw / (1 - γ)

    params.cache.γ_SEE .= γ
    params.cache.νew[1] = νew
    params.cache.νew[2] = νew
    params.cache.νew[3] = 0.0
    params.cache.νew[4] = 0.0

    Iiw = HallThruster.wall_ion_current(sheath_model, U, params, 2, 1)
    Iew = HallThruster.wall_electron_current(sheath_model, U, params, 2)
    @test Iew ≈ Iiw / (1 - γ)
    @test Iiw ≈ νiw * HallThruster.e * V_cell * ne
    @test Iew ≈ νew * HallThruster.e * V_cell * ne

    @test HallThruster.freq_electron_wall(sheath_model, U, params, 2) ≈ νew
    @test HallThruster.freq_electron_wall(sheath_model, U, params, 3) ≈ 0.0

    @test HallThruster.wall_power_loss(sheath_model, U, params, 2) ≈ νew * (2 * (1 - 0.5 * BN.σ₀) * Tev + (1 - γ) * Vs)/γ
    @test HallThruster.wall_power_loss(sheath_model, U, params, 4) ≈ 0.0
end

@testset "Ion wall losses" begin
    config = (;thruster = HallThruster.SPT_100, propellant = HallThruster.Krypton, ncharge = 2, transition_function = HallThruster.StepFunction())
    geom = HallThruster.SPT_100.geometry
    Δr = geom.outer_radius - geom.inner_radius

    Tn = 300.0
    Tev = 3.0

    α = 0.8
    mi = config.propellant.m

    γ_SEE_max = 1 - 8.3 * sqrt(HallThruster.me/mi)

    γ = HallThruster.SEE_yield(HallThruster.BoronNitride, Tev, γ_SEE_max)
    νiw = α * sqrt(HallThruster.e * Tev / mi) / Δr
    νew = νiw * γ / (1 - γ)

    L_ch = HallThruster.geometry_SPT_100.channel_length
    A_ch = HallThruster.channel_area(HallThruster.SPT_100)

    grid = HallThruster.generate_grid(HallThruster.geometry_SPT_100, 2, (0, 2 * L_ch))

    z_cell = grid.cell_centers
    z_edge = grid.edges

    # Fill up cell lengths and magnetic field vectors
    Δz_cell = zeros(length(z_cell))
    Δz_edge = zeros(length(z_edge))
    for (i, z) in enumerate(grid.cell_centers)
        if firstindex(z_cell) < i < lastindex(z_cell)
            Δz_cell[i] = z_edge[HallThruster.right_edge(i)] - z_edge[HallThruster.left_edge(i)]
        elseif i == firstindex(z_cell)
            Δz_cell[i] = z_edge[begin+1] - z_edge[begin]
        elseif i == lastindex(z_cell)
            Δz_cell[i] = z_edge[end] - z_edge[end-1]
        end
    end

    for i in eachindex(z_edge)
        Δz_edge[i] = z_cell[i+1] - z_cell[i]
    end

    fluids = [
        HallThruster.Fluid(Xenon(0), HallThruster.ContinuityOnly(u = 300.0, T = Tn)),
        HallThruster.Fluid(Xenon(1), HallThruster.IsothermalEuler(T = Tn)),
        HallThruster.Fluid(Xenon(2), HallThruster.IsothermalEuler(T = Tn)),
    ]

    Tev = 4.0
    ne = 1.3e18
    ni_1 = 7e17
    ni_2 = 3e17

    cache = (;
        ne = [ne, ne, ne, ne], Tev = [Tev, Tev, Tev, Tev],
        Z_eff = [1.0, 1.0, 1.0, 1.0], ni = [ni_1 ni_1 ni_1 ni_1; ni_2 ni_2 ni_2 ni_2],
        γ_SEE = [0.0, 0.0, 0.0, 0.0],
        νew = [νew, νew, 0.0, 0.0]
    )

    index = (ρn = 1, ρi = [2, 4], ρiui = [3, 5], nϵ = 6)

    u_bohm = sqrt(HallThruster.e * Tev / mi)
    u_thermal_n = sqrt(HallThruster.kB * Tn / 2 / π / mi)

    ρn = 10 * ne * mi

    ρi_1 = ni_1 * mi
    ρi_2 = ni_2 * mi
    ui_1 = 10.0
    ui_2 = sqrt(2) * 10.0

    nϵ = ne * Tev * 3/2

    U = [
        ρn          ρn              ρn          ρn
        ρi_1        ρi_1            ρi_1        ρi_1
        ρi_1 * ui_1 ρi_1 * ui_1     ρi_1 * ui_1 ρi_1 * ui_1
        ρi_2        ρi_2            ρi_2        ρi_2
        ρi_2 * ui_2 ρi_2 * ui_2     ρi_2 * ui_2 ρi_2 * ui_2
        nϵ          nϵ              nϵ          nϵ
    ]

    dU = zeros(size(U))

    z_edge = grid.edges
    z_cell = grid.cell_centers
    γ_SEE_max = 1 - 8.3 * sqrt(HallThruster.me/mi)
    base_params = (;cache, L_ch, A_ch, fluids, z_cell, z_edge, index, Δz_cell, Δz_edge, γ_SEE_max)

    config_no_losses = (;config..., wall_loss_model = HallThruster.NoWallLosses())
    params_no_losses = (;base_params..., config = config_no_losses)

    for i in 1:4
        cache.Z_eff[i] = HallThruster.compute_Z_eff(U, params_no_losses, i)
    end

    @test cache.Z_eff[1] ≈ (ni_1 + 2 * ni_2) / (ni_1 + ni_2)

    u_bohm_1 = u_bohm
    u_bohm_2 = sqrt(2) * u_bohm

    # Test 1: no wall losses
    HallThruster.apply_ion_wall_losses!(dU, U, params_no_losses, 2)
    HallThruster.apply_ion_wall_losses!(dU, U, params_no_losses, 3)

    @test all(dU .≈ 0.0)

    # Test 2: LANDMARK wall losses

    αin, αout = 1.0, 1.0
    constant_sheath = HallThruster.ConstantSheathPotential(-20, αin, αout)
    config_constant_sheath = (;config..., wall_loss_model = constant_sheath)
    params_constant_sheath = (;base_params..., config = config_constant_sheath)

    i = 2
    Iiw_1 = HallThruster.wall_ion_current(constant_sheath, U, params_constant_sheath, i, 1)
    Iiw_2 = HallThruster.wall_ion_current(constant_sheath, U, params_constant_sheath, i, 2)
    Iew = HallThruster.wall_electron_current(constant_sheath, U, params_constant_sheath, i)

    # Ion and electron wall currents are equivalent inside of channel
    @test Iiw_1 + Iiw_2 == Iew

    dU .= 0.0
    # Check that wall losses work correctly
    HallThruster.apply_ion_wall_losses!(dU, U, params_constant_sheath, i)

    # Neutrals should recombine at walls
    @test dU[index.ρn[1], i] ≈ -(dU[index.ρi[1], i] + dU[index.ρi[2], i])

    # Rate of ion loss is equal to ni νiw
    Δz = z_edge[2] - z_edge[1]
    V_cell = A_ch * Δz
    u_bohm = sqrt(HallThruster.e * Tev / mi)
    νiw =  u_bohm / Δr

    @test dU[index.ρi[1], i] ≈ -U[index.ρi[1], i] * νiw
    @test dU[index.ρi[2], i] ≈ -U[index.ρi[2], i] * νiw * sqrt(2)

    # ion momentum loss is equal to Iiw / e / V_cell * ui
    @test dU[index.ρiui[1], i] ≈ -U[index.ρiui[1], i] * νiw
    @test dU[index.ρiui[2], i] ≈ -U[index.ρiui[2], i] * νiw * sqrt(2)

    # No wall current outside of channel
    dU .= 0.0
    i = 3
    @test HallThruster.wall_ion_current(constant_sheath, U, params_constant_sheath, i, 1) == 0.0
    @test HallThruster.wall_ion_current(constant_sheath, U, params_constant_sheath, i, 2) == 0.0
    @test HallThruster.wall_electron_current(constant_sheath, U, params_constant_sheath, i) == 0.0

    HallThruster.apply_ion_wall_losses!(dU, U, params_constant_sheath, 3)
    @test all(dU[:, 3] .≈ 0.0)

    # Test 3: Self-consistent wall sheath
    dU .= 0.0
    material = HallThruster.BNSiO2
    wall_sheath = HallThruster.WallSheath(material, α)
    config_wall_sheath = (;config..., wall_loss_model = wall_sheath)
    params_wall_sheath = (;base_params..., config = config_wall_sheath)

    γ = HallThruster.SEE_yield(material, Tev, γ_SEE_max)
    νiw = α * sqrt(HallThruster.e * Tev / mi) / Δr
    νew = νiw * γ / (1 - γ)

    params_wall_sheath.cache.νew[1:2] .= νew
    i = 2
    params_wall_sheath.cache.γ_SEE[i] = γ
    Iiw_1 = HallThruster.wall_ion_current(wall_sheath, U, params_wall_sheath, i, 1)
    Iiw_2 = HallThruster.wall_ion_current(wall_sheath, U, params_wall_sheath, i, 2)
    Iew = HallThruster.wall_electron_current(wall_sheath, U, params_wall_sheath, i)

    @test Iiw_1 ≈ Iew * ni_1 / ne * (1 - γ)
    @test Iiw_2 ≈ 2 * Iew * ni_2 / ne * (1 - γ)
    @test Iew ≈ inv(1 - γ) * (Iiw_1 + Iiw_2)

    HallThruster.apply_ion_wall_losses!(dU, U, params_wall_sheath, i)

    # Neutrals should recombine at walls
    @test dU[index.ρn[1], i] ≈ -(dU[index.ρi[1], i] + dU[index.ρi[2], i])

    # Rate of ion loss is equal to Iiw / e / V_cell
    @test dU[index.ρi[1], i] ≈ -νiw * U[index.ρi[1], i]
    @test dU[index.ρi[2], i] ≈ -νiw * U[index.ρi[2], i] * sqrt(2)

    # ion momentum loss is equal to Iiw / e / V_cell * ui
    @test dU[index.ρiui[1], i] ≈ -νiw * U[index.ρiui[1], i]
    @test dU[index.ρiui[2], i] ≈ -νiw * U[index.ρiui[2], i] * sqrt(2)

    # No wall losses in plume
    i = 3
    @test HallThruster.wall_ion_current(wall_sheath, U, params_wall_sheath, i, 1) == 0.0
    @test HallThruster.wall_ion_current(wall_sheath, U, params_wall_sheath, i, 2) == 0.0
    @test HallThruster.wall_electron_current(wall_sheath, U, params_wall_sheath, i) == 0.0

    HallThruster.apply_ion_wall_losses!(dU, U, params_wall_sheath, 3)
    @test all(dU[:, 3] .≈ 0.0)
end
