@testset "Ionization" begin
    Xe_0    = HallThruster.Xenon(0)
    Xe_I    = HallThruster.Xenon(1)
    Xe_II   = HallThruster.Xenon(2)
    Xe_III  = HallThruster.Xenon(3)
    Xe_IV   = HallThruster.Xenon(4)

    rxn_0_I = HallThruster.IonizationReaction(Xe_0, Xe_I, Te -> 0.0)
    rxn_0_II = HallThruster.IonizationReaction(Xe_0, Xe_II, Te -> 0.0)
    rxn_0_III = HallThruster.IonizationReaction(Xe_0, Xe_III, Te -> 0.0)
    rxn_I_III = HallThruster.IonizationReaction(Xe_I, Xe_III, Te -> 0.0)
    @test repr(rxn_0_I) == "e- + Xe -> 2e- + Xe+"
    @test repr(rxn_0_II) == "e- + Xe -> 3e- + Xe2+"
    @test repr(rxn_0_III) == "e- + Xe -> 4e- + Xe3+"
    @test repr(rxn_I_III) == "e- + Xe+ -> 3e- + Xe3+"

    @test HallThruster.rate_coeff_filename(Xe_0, Xe_II, "ionization") == joinpath(HallThruster.REACTION_FOLDER, "ionization_Xe_Xe2+.dat")

    @test_throws(ArgumentError, HallThruster.load_ionization_reaction(Xe_II, Xe_0))
    @test !isnothing(HallThruster.load_ionization_reaction(Xe_0, Xe_II))

    struct MyModel <: HallThruster.IonizationModel end

    mymodel = MyModel()

    @test HallThruster.supported_gases(mymodel) == HallThruster.Gas[]
    @test HallThruster.maximum_charge_state(mymodel) == 0

    @test_throws ArgumentError HallThruster.load_ionization_reactions(mymodel, [Xe_0])

    landmark_lut = HallThruster.LandmarkIonizationLUT()
    bolsig_lut = HallThruster.BolsigIonizationLUT()
    bolsig_fit = HallThruster.BolsigIonizationFit()

    @test HallThruster.supported_gases(landmark_lut) == [HallThruster.Xenon]
    @test HallThruster.maximum_charge_state(landmark_lut) == 1

    @test HallThruster.supported_gases(bolsig_lut) == [HallThruster.Xenon, HallThruster.Krypton]
    @test HallThruster.maximum_charge_state(bolsig_lut) == 3

    @test HallThruster.supported_gases(bolsig_fit) == [HallThruster.Xenon]
    @test HallThruster.maximum_charge_state(bolsig_fit) == 3

    Ar_0 = HallThruster.Argon(0)
    Ar_I = HallThruster.Argon(1)

    @test_throws ArgumentError HallThruster._load_ionization_reactions(bolsig_fit, [Ar_0, Ar_I])
    @test_throws ArgumentError HallThruster._load_ionization_reactions(bolsig_lut, [Ar_0, Ar_I])
    @test_throws ArgumentError HallThruster._load_ionization_reactions(landmark_lut, [Ar_0, Ar_I])

    @test_throws ArgumentError HallThruster._load_ionization_reactions(landmark_lut, [Xe_0, Xe_I, Xe_II])
    @test_throws ArgumentError HallThruster._load_ionization_reactions(landmark_lut, [Xe_0, Xe_I, Xe_II, Xe_III])
    @test_throws ArgumentError HallThruster._load_ionization_reactions(bolsig_fit, [Xe_0, Xe_I, Xe_II, Xe_III, Xe_IV])
    @test_throws ArgumentError HallThruster._load_ionization_reactions(bolsig_lut, [Xe_0, Xe_I, Xe_II, Xe_III, Xe_IV])

    landmark_rxns = HallThruster._load_ionization_reactions(landmark_lut, [Xe_0, Xe_I])
    bolsig_rxns = HallThruster._load_ionization_reactions(bolsig_lut, [Xe_0, Xe_I, Xe_II, Xe_III])
    bolsig_fit_rxns = HallThruster._load_ionization_reactions(bolsig_fit, [Xe_0, Xe_I, Xe_II, Xe_III])

    @test length(landmark_rxns) == 1
    @test length(bolsig_rxns) == 6
    @test length(bolsig_fit_rxns) == 6

    @test landmark_rxns[1].rate_coeff(19) ≈ 5.6880E-14

    @test count(rxn -> rxn.reactant == Xe_0, bolsig_rxns) == 3
    @test count(rxn -> rxn.reactant == Xe_I, bolsig_rxns) == 2
    @test count(rxn -> rxn.reactant == Xe_II, bolsig_rxns) == 1
    @test count(rxn -> rxn.reactant == Xe_III, bolsig_rxns) == 0
    @test count(rxn -> rxn.product == Xe_0, bolsig_rxns) == 0
    @test count(rxn -> rxn.product == Xe_I, bolsig_rxns) == 1
    @test count(rxn -> rxn.product == Xe_II, bolsig_rxns) == 2
    @test count(rxn -> rxn.product == Xe_III, bolsig_rxns) == 3

    @test_throws(ArgumentError, HallThruster.ionization_fits_Xe(0))
    @test_throws(ArgumentError, HallThruster.ionization_fits_Xe(4))
    @test length(HallThruster.ionization_fits_Xe(1)) == 1
    @test length(HallThruster.ionization_fits_Xe(2)) == 3
    @test length(HallThruster.ionization_fits_Xe(3)) == 6
end