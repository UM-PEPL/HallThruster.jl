Xe_0 = HallThruster.Species(HallThruster.Xenon, 0)
Xe_I = HallThruster.Species(HallThruster.Xenon, 1)
Xe_II = HallThruster.Species(HallThruster.Xenon, 2)
Xe_III = HallThruster.Species(HallThruster.Xenon, 3)

rxn_0_I = HallThruster.IonizationReaction(Xe_0, Xe_I, Te -> 0.0)
rxn_0_II = HallThruster.IonizationReaction(Xe_0, Xe_II, Te -> 0.0)
rxn_0_III = HallThruster.IonizationReaction(Xe_0, Xe_III, Te -> 0.0)
rxn_I_III = HallThruster.IonizationReaction(Xe_I, Xe_III, Te -> 0.0)
@test repr(rxn_0_I) == "e- + Xe -> 2e- + Xe+"
@test repr(rxn_0_II) == "e- + Xe -> 3e- + Xe2+"
@test repr(rxn_0_III) == "e- + Xe -> 4e- + Xe3+"
@test repr(rxn_I_III) == "e- + Xe+ -> 3e- + Xe3+"

@test HallThruster.rate_coeff_filename(Xe_0, Xe_II, "ionization") == "ionization_Xe_Xe2+.dat"

@test_throws(ArgumentError, HallThruster.load_ionization_reaction(Xe_II, Xe_0))
@test !isnothing(HallThruster.load_ionization_reaction(Xe_0, Xe_II))