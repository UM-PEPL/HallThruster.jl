Tev = 30 #[eV]
m = HallThruster.Xenon.m #
Te = Tev*HallThruster.e/HallThruster.kB #
ne = 1e18 #[#/m^3]
nn = 0.5e18 #[#/m^3]
B = 1.0
ν_an = 0.0
σ_en = 6.6e-19*((Tev/4 - 0.1)/(1 + (Tev/4)^1.6)) #[m^2]
@test σ_en ≈ HallThruster.σ_en(Tev)
ln_λ = 24 - 0.5*log(1e-6*ne/Tev^2)
@test ln_λ ≈ HallThruster.coulomb_logarithm(ne, Tev)
Tev = 9
ln_λ = 23 - 0.5*log(1e-6*ne/Tev^3)
@test ln_λ ≈ HallThruster.coulomb_logarithm(ne, Tev)
ν_c = σ_en*nn*sqrt(8*HallThruster.kB*Te/pi/m) + 2.9e-12*ne*ln_λ/(Tev)^1.5
#@test ν_c ≈ HallThruster.get_v_c(Tev, nn, ne, m) #can't pass if Landmark set
μ_e = HallThruster.e/(HallThruster.mₑ * ν_c)/(1+(HallThruster.e*B/(HallThruster.mₑ*ν_c))^2)
#@test μ_e ≈ HallThruster.cf_electron_transport(ν_an, ν_c, B) can't pass if Landmark set