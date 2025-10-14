using HallThruster: HallThruster as het
using DelimitedFiles

Xe_0 = het.Xenon(0)
Xe_I = het.Xenon(1)
Xe_II = het.Xenon(2)
Xe_III = het.Xenon(3)
Xe_IV = het.Xenon(4)

r = zeros(256)
rxn_0_I = het.ElectronImpactReaction(0.0, Xe_0, [Xe_I], r)
rxn_0_II = het.ElectronImpactReaction(0.0, Xe_0, [Xe_II], r)
rxn_0_III = het.ElectronImpactReaction(0.0, Xe_0, [Xe_III], r)
rxn_I_III = het.ElectronImpactReaction(0.0, Xe_I, [Xe_III], r)
@test repr(rxn_0_I) == "e(-) + Xe -> 2e(-) + Xe(+)"
@test repr(rxn_0_II) == "e(-) + Xe -> 3e(-) + Xe(2+)"
@test repr(rxn_0_III) == "e(-) + Xe -> 4e(-) + Xe(3+)"
@test repr(rxn_I_III) == "e(-) + Xe(+) -> 3e(-) + Xe(3+)"

@test het.rate_coeff_filename(Xe_0, Xe_II, "ionization") ==
    joinpath(het.REACTION_FOLDER, "ionization_Xe_Xe2+.dat")
@test het.rate_coeff_filename(Xe_0, Xe_II, "excitation") ==
    joinpath(het.REACTION_FOLDER, "excitation_Xe_Xe2+.dat")
@test het.rate_coeff_filename(Xe_0, nothing, "excitation") ==
    joinpath(het.REACTION_FOLDER, "excitation_Xe.dat")

Bi_0 = het.Bismuth(0)
Bi_I = het.Bismuth(1)

# Test behavior of landmark lookup
@test_throws ArgumentError het.load_electron_impact_reactions(:Landmark, [Bi_0, Bi_I])
@test_throws ArgumentError het.load_electron_impact_reactions(
    :Landmark, [Xe_0, Xe_I, Xe_II],
)
@test_throws ArgumentError het.load_electron_impact_reactions(
    :Landmark, [Xe_0, Xe_I, Xe_II, Xe_III],
)

landmark_rxns = het.load_electron_impact_reactions(:Landmark, [Xe_0, Xe_I])
@test length(landmark_rxns) == 1
@test het.rate_coeff(landmark_rxns[1], 19.0) ≈ 5.69e-14

# Test behavior of general lookup
@test_throws ArgumentError het.load_electron_impact_reactions(:Lookup, [Bi_0, Bi_I])
@test_throws ArgumentError het.load_electron_impact_reactions(
    :Lookup, [Xe_0, Xe_I, Xe_II, Xe_III, Xe_IV],
)
lookup_rxns = het.load_electron_impact_reactions(:Lookup, [Xe_0, Xe_I, Xe_II, Xe_III])
@test length(lookup_rxns) == 6
@test count(rxn -> rxn.reactant == Xe_0, lookup_rxns) == 3
@test count(rxn -> rxn.reactant == Xe_I, lookup_rxns) == 2
@test count(rxn -> rxn.reactant == Xe_II, lookup_rxns) == 1
@test count(rxn -> rxn.reactant == Xe_III, lookup_rxns) == 0
@test count(rxn -> rxn.products[] == Xe_I, lookup_rxns) == 1
@test count(rxn -> rxn.products[] == Xe_0, lookup_rxns) == 0
@test count(rxn -> rxn.products[] == Xe_II, lookup_rxns) == 2
@test count(rxn -> rxn.products[] == Xe_III, lookup_rxns) == 3

# Test behavior of user-provided ionization reactions
directories = [
    joinpath(
        het.PACKAGE_ROOT, "test", "unit_tests", "reaction_tests",
    ),
]

lookup_2_rxns = het.load_electron_impact_reactions(:Lookup, [Bi_0, Bi_I]; directories)
@test length(lookup_2_rxns) == 1
@test het.rate_coeff(lookup_2_rxns[1], 0.3878e-1) |> abs <
    eps(Float64)
@test lookup_2_rxns[1].energy == 13.0
lookup_2_rxns_Xe = het.load_electron_impact_reactions(
    :Lookup, [Xe_0, Xe_I, Xe_II, Xe_III]; directories,
)
@test length(lookup_2_rxns_Xe) == 6
@test het.rate_coeff(lookup_2_rxns_Xe[1], 1.0) |> abs <
    eps(Float64)
@test het.rate_coeff(lookup_2_rxns_Xe[1], 1.0) !=
    het.rate_coeff(lookup_rxns[1], 1.0)
@test lookup_2_rxns_Xe[1].energy == 15.0

# Excitation reactions
ex_1 = het.ExcitationReaction(0.0, Xe_0, r)
@test repr(ex_1) == "e(-) + Xe -> e(-) + Xe(*)"

@test isempty(het.load_excitation_reactions(:Lookup, [Bi_0]))
ex_rxns = het.load_excitation_reactions(:Lookup, [Xe_0])
@test length(ex_rxns) == 1
@test ex_rxns[1].energy == 8.32
@test het.rate_coeff(ex_rxns[1], 10.0) ≈ 2.358251483e-14

@test !isempty(het.load_excitation_reactions(:Lookup, [Bi_0]; directories))
ex_rxns = het.load_excitation_reactions(:Lookup, [Bi_0]; directories)
@test ex_rxns[1].energy == 10.0
@test het.rate_coeff(ex_rxns[1], 1.0) |> abs < eps(Float64)

ex_landmark_rxn = het.load_excitation_reactions(:Landmark, [Xe_0])[1]
@test ex_landmark_rxn.energy == 8.32

landmark_table = readdlm(
    joinpath(het.PACKAGE_ROOT, "landmark", "landmark_rates.csv"),
    ',', skipstart = 1,
)
loss_itp = het.LinearInterpolation(landmark_table[:, 1], landmark_table[:, 3])
iz_landmark_rxn = landmark_rxns[1]

@test 8.32 * het.rate_coeff(ex_landmark_rxn, 14.32) +
    12.12 * het.rate_coeff(iz_landmark_rxn, 14.32) ≈
    loss_itp(14.32)

# More complex reactions
rxn_test = het.ElectronImpactReaction(0.0, het.MolecularNitrogen(0), [het.Nitrogen(0), het.Nitrogen(0)], r)
@test repr(rxn_test) == "e(-) + N2 -> e(-) + N + N"
