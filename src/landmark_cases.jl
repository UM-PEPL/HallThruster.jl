"""
    SPT_100
Geometry of the SPT_100 thruster
"""
const SPT_100 = (domain=(0.0, 0.05), channel_length=0.025, inner_radius=0.0345,
                 outer_radius=0.05)

"""
    B_field(B_max::Float64, z::Float64, L_ch::Float64)

defines magnetic field as a function of position. 
"""
function B_field(B_max, z, L_ch) #same in Landmark and in FFM model Hara
    B = if z < L_ch
        B_max * exp(-0.5 * ((z - L_ch) / (0.011))^2) #for SPT_100
    else
        B_max * exp(-0.5 * ((z - L_ch) / (0.018))^2)
    end
    return B
end


const LANDMARK_MAGNETIC_FIELD = let z = 0:0.001:0.05
    MagneticField(collect(z), B_field.(0.015, z, 0.025))
end

function precompute_bfield!(B, zs)
    B_max = 0.015
    L_ch = 0.025
    for (i, z) in enumerate(zs)
        B[i] = B_field(B_max, z, L_ch)
    end
end

function landmark_case(case_num; implicit_energy = false, dt = 1e-9, nsave=100, sim_time = 5e-4, ncells=50, adaptive=true)
    νϵ = if case_num == 1
        1e7
    elseif case_num == 2
        0.5e7
    elseif case_num == 3
        0.4e7
    end

    landmark = HallThrusterConfig(
        name = "Landmark $(case_num)",
        propellant = Xenon,
        ncells = ncells,
        ncharge = 1,
        magnetic_field = LANDMARK_MAGNETIC_FIELD,
        dt = dt,
        dtmax = dt,
        nsave = nsave,
        simulation_time = sim_time,
        implicit_energy = implicit_energy,
        cathode_Te = 3.0,
        anode_Te = 3.0,
        cathode_potential = 0.0,
        anode_potential = 300.0,
        geometry = SPT_100,
        anom_model = TwoZoneBohm(),
        radial_loss_coefficients = (νϵ, 1e7),
        wall_collision_frequencies = (1e7, 1e7),
        energy_equation = :LANDMARK,
        ionization_coeffs = :LANDMARK,
        ion_temperature = 500.0,
        ion_diffusion_coeff = 0.0,#0.5e-3,
        coupled_method = true,
        neutral_velocity = 150.0,
        neutral_temperature = 300.0,
        adaptive = adaptive,
        anode_mass_flow_rate = 5e-6,
    )

    @time sol = run_simulation(landmark)

    return sol
end