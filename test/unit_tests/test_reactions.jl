@testset "Ionization" begin
    Xe_0    = HallThruster.Xenon(0)
    Xe_I    = HallThruster.Xenon(1)
    Xe_II   = HallThruster.Xenon(2)
    Xe_III  = HallThruster.Xenon(3)
    Xe_IV   = HallThruster.Xenon(4)

    rxn_0_I = HallThruster.IonizationReaction(0.0, Xe_0, Xe_I, Te -> 0.0)
    rxn_0_II = HallThruster.IonizationReaction(0.0, Xe_0, Xe_II, Te -> 0.0)
    rxn_0_III = HallThruster.IonizationReaction(0.0, Xe_0, Xe_III, Te -> 0.0)
    rxn_I_III = HallThruster.IonizationReaction(0.0, Xe_I, Xe_III, Te -> 0.0)
    @test repr(rxn_0_I) == "e- + Xe -> 2e- + Xe+"
    @test repr(rxn_0_II) == "e- + Xe -> 3e- + Xe2+"
    @test repr(rxn_0_III) == "e- + Xe -> 4e- + Xe3+"
    @test repr(rxn_I_III) == "e- + Xe+ -> 3e- + Xe3+"

    @test HallThruster.rate_coeff_filename(Xe_0, Xe_II, "ionization") == joinpath(HallThruster.REACTION_FOLDER, "ionization_Xe_Xe2+.dat")
    @test HallThruster.rate_coeff_filename(Xe_0, Xe_II, "excitation") == joinpath(HallThruster.REACTION_FOLDER, "excitation_Xe_Xe2+.dat")
    @test HallThruster.rate_coeff_filename(Xe_0, nothing, "excitation") == joinpath(HallThruster.REACTION_FOLDER, "excitation_Xe.dat")

    @test_throws(ArgumentError, HallThruster.load_rate_coeffs(Xe_II, Xe_0, "ionization"))
    @test !isnothing(HallThruster.load_rate_coeffs(Xe_0, Xe_II, "ionization"))

    struct MyModel <: HallThruster.IonizationModel end

    mymodel = MyModel()

    @test HallThruster.supported_gases(mymodel) == HallThruster.Gas[]
    @test HallThruster.maximum_charge_state(mymodel) == 0

    @test_throws ArgumentError HallThruster.load_reactions(mymodel, [Xe_0])

    landmark_lut = HallThruster.LandmarkIonizationLookup()
    lookup = HallThruster.IonizationLookup()

    @test HallThruster.supported_gases(landmark_lut) == [HallThruster.Xenon]
    @test HallThruster.maximum_charge_state(landmark_lut) == 1

    @test HallThruster.supported_gases(lookup) == HallThruster.Gas[]
    @test HallThruster.maximum_charge_state(lookup) == 0

    Ar_0 = HallThruster.Argon(0)
    Ar_I = HallThruster.Argon(1)

    # Test behavior of landmark lookup
    @test_throws ArgumentError HallThruster._load_reactions(landmark_lut, [Ar_0, Ar_I])
    @test_throws ArgumentError HallThruster._load_reactions(landmark_lut, [Xe_0, Xe_I, Xe_II])
    @test_throws ArgumentError HallThruster._load_reactions(landmark_lut, [Xe_0, Xe_I, Xe_II, Xe_III])
    landmark_rxns = HallThruster._load_reactions(landmark_lut, [Xe_0, Xe_I])
    @test length(landmark_rxns) == 1
    @test landmark_rxns[1].rate_coeff(19) ≈ 5.690E-14

    # Test behavior of general lookup
    @test_throws ArgumentError HallThruster.load_reactions(lookup, [Ar_0, Ar_I])
    @test_throws ArgumentError HallThruster.load_reactions(lookup, [Xe_0, Xe_I, Xe_II, Xe_III, Xe_IV])
    lookup_rxns = HallThruster._load_reactions(lookup, [Xe_0, Xe_I, Xe_II, Xe_III])
    @test length(lookup_rxns) == 6
    @test count(rxn -> rxn.reactant == Xe_0, lookup_rxns) == 3
    @test count(rxn -> rxn.reactant == Xe_I, lookup_rxns) == 2
    @test count(rxn -> rxn.reactant == Xe_II, lookup_rxns) == 1
    @test count(rxn -> rxn.reactant == Xe_III, lookup_rxns) == 0
    @test count(rxn -> rxn.product == Xe_I, lookup_rxns) == 1
    @test count(rxn -> rxn.product == Xe_0, lookup_rxns) == 0
    @test count(rxn -> rxn.product == Xe_II, lookup_rxns) == 2
    @test count(rxn -> rxn.product == Xe_III, lookup_rxns) == 3

    # Test behavior of user-provided ionization reactions
    lookup_2 = HallThruster.IonizationLookup([joinpath(HallThruster.PACKAGE_ROOT, "test", "unit_tests", "reaction_tests")])
    lookup_2_rxns = HallThruster._load_reactions(lookup_2, [Ar_0, Ar_I])
    @test length(lookup_2_rxns) == 1
    @test lookup_2_rxns[1].rate_coeff(0.3878e-01) ≈ 0.0
    @test lookup_2_rxns[1].energy == 13.0
    lookup_2_rxns_Xe = HallThruster._load_reactions(lookup_2, [Xe_0, Xe_I, Xe_II, Xe_III])
    @test length(lookup_2_rxns_Xe) == 6
    @test lookup_2_rxns_Xe[1].rate_coeff(1.0) ≈ 0.0
    @test lookup_2_rxns_Xe[1].rate_coeff(1.0) != lookup_rxns[1].rate_coeff(1.0)
    @test lookup_2_rxns_Xe[1].energy == 15.0

    # Excitation reactions
    ex_1 = HallThruster.ExcitationReaction(0.0, Xe_0, Te -> 0.0)
    @test repr(ex_1) == "e- + Xe -> e- + Xe*"

    ex_lookup = HallThruster.ExcitationLookup()
    @test isempty(HallThruster._load_reactions(ex_lookup, [Ar_0]))
    ex_rxns = HallThruster._load_reactions(ex_lookup, [Xe_0])
    @test length(ex_rxns) == 1
    @test ex_rxns[1].energy == 8.32
    @test ex_rxns[1].rate_coeff(10) ≈ 2.3582514835893516e-14

    ex_lookup_2 = HallThruster.ExcitationLookup([joinpath(HallThruster.PACKAGE_ROOT, "test", "unit_tests", "reaction_tests")])
    @test !isempty(HallThruster._load_reactions(ex_lookup_2, [Ar_0]))
    ex_rxns = HallThruster._load_reactions(ex_lookup_2, [Ar_0])
    @test ex_rxns[1].energy == 10.0
    @test ex_rxns[1].rate_coeff(1.0) ≈ 0.0

    ex_lookup_landmark = HallThruster.LandmarkExcitationLookup()
    ex_landmark_rxn = HallThruster.load_reactions(ex_lookup_landmark, [Xe_0])[1]
    @test ex_landmark_rxn.energy == 8.32

    landmark_table = readdlm(joinpath(HallThruster.PACKAGE_ROOT, "landmark", "landmark_rates.csv"), ',', skipstart = 1)
    loss_itp = HallThruster.LinearInterpolation(landmark_table[:, 1], landmark_table[:, 3])
    iz_landmark_rxn = landmark_rxns[1]

    @test 8.32 * ex_landmark_rxn.rate_coeff(14.32) + 12.12 * iz_landmark_rxn.rate_coeff(14.32) ≈ loss_itp(14.32)
end