@testset "Ionization" begin
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

    let landmark = HallThruster.load_landmark()
        @test landmark.rate_coeff(19) ≈ 5.6880E-14
        @test landmark.loss_coeff(19) ≈ 1.1822E-12
    end

    let species = [Xe_0, Xe_I, Xe_II, Xe_III]
        reactions = HallThruster.load_ionization_reactions(species)

        @test length(reactions) == 6
        @test count(rxn -> rxn.reactant == Xe_0, reactions) == 3
        @test count(rxn -> rxn.reactant == Xe_I, reactions) == 2
        @test count(rxn -> rxn.reactant == Xe_II, reactions) == 1
        @test count(rxn -> rxn.reactant == Xe_III, reactions) == 0
        @test count(rxn -> rxn.product == Xe_0, reactions) == 0
        @test count(rxn -> rxn.product == Xe_I, reactions) == 1
        @test count(rxn -> rxn.product == Xe_II, reactions) == 2
        @test count(rxn -> rxn.product == Xe_III, reactions) == 3
    end

    @test_throws(ArgumentError, HallThruster.ionization_fits_Xe(0))
    @test_throws(ArgumentError, HallThruster.ionization_fits_Xe(4))
    @test length(HallThruster.ionization_fits_Xe(1)) == 1
    @test length(HallThruster.ionization_fits_Xe(2)) == 3
    @test length(HallThruster.ionization_fits_Xe(3)) == 6
end