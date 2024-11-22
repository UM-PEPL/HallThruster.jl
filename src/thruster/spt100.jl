const spt100_geometry = Geometry1D(;
    inner_radius = 0.0345, outer_radius = 0.05, channel_length = 0.025,
)

function spt100_analytic_field()
    L_ch = spt100_geometry.channel_length
    B_max = 0.015 # T
    zs = LinRange(0, 4 * L_ch, 256)
    Bs = zeros(length(zs))
    for (i, z) in enumerate(zs)
        Bs[i] = if z < L_ch
            B_max * exp(-0.5 * ((z - L_ch) / (0.011))^2)
        else
            B_max * exp(-0.5 * ((z - L_ch) / (0.018))^2)
        end
    end
    return MagneticField("SPT-100 analytic field", zs, Bs)
end

const SPT_100 = Thruster(;
    name = "SPT-100",
    geometry = spt100_geometry,
    magnetic_field = spt100_analytic_field(),
    shielded = false,
)
