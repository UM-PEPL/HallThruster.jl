@inline function σ_en(Tev) #Te in eV, from intro to EP, 3.6-13
    return 6.6e-19*((Tev/4-0.1)/(1+(Tev/4)^1.6)) #[m^2] 
end

@inline function ln_λ(ne, Tev) #from intro to EP, 3.6-15, or just assume constant 15-25
    return 23 - 0.5*log(ne*1e-6/Tev^3)
end

@inline function get_v_c(Tev, nn, ne, m) #classical momentum transfer collision frequency, v_ei + v_en, nn neutral number density, ne electron number density, m ion mass
    v_en = σ_en(Tev)*nn*sqrt(8*e*Tev/pi/m) # intro to EP, 3.6-12
    v_ei = 2.9e-12*ne*ln_λ(ne, Tev)/Tev^1.5 #intro to EP, 3.6-14
    return v_en + v_ei #2.5e-13*nn #from Hara paper, is similar to formula for v_ei from intro to EP
    #return 2.5e-13*nn
end

@inline function get_v_an(z, B, L_ch) #anomalous momentum transfer collision frequency
    ωce = e * B / mₑ
    #νan = ωce/16
    νan = if z < L_ch
        ωce/160
    else
        ωce/16
    end
    return νan
end

function B_field(B_max, z, L_ch)
    B = if z < L_ch
        B_max*exp(-0.5*((z-L_ch)/(0.011))^2) #for SPT_100
    else
        B_max*exp(-0.5*((z-L_ch)/(0.018))^2)
    end
    return B
end

function Te_func(z, L_ch)
    return 30 * exp(-(2 * (z - L_ch) / 0.033)^2)
end

@inline function cf_electron_transport(v_an, v_c, B)
    vₑ = v_an + v_c
    Ω = e*B/(mₑ*vₑ)
    return e/(mₑ*vₑ*(1+Ω^2))
end 

function electron_pressure(ne, Tev)
    return e * ne * Tev
end