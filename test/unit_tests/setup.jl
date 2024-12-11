using Unitful

function test_setup()
	@testset "Simulation setup" begin
		geom = het.Geometry1D(
			channel_length = 2.5u"cm",
			inner_radius = 3.45u"cm",
			outer_radius = 5u"cm",
		)

		@test geom.channel_length ≈ 0.025
		@test geom.inner_radius ≈ 0.0345
		@test geom.outer_radius ≈ 0.05

		simparams = het.SimParams(
			grid = het.EvenGrid(100),
			dt = 5u"ns",
			duration = 1u"ms",
			num_save = 1000,
			verbose = false,
		)
		@test simparams.dt ≈ 5e-9
		@test simparams.duration ≈ 1e-3

		config = het.Config(
			thruster = het.SPT_100,
			domain = (0.0u"cm", 8.0u"cm"),
			discharge_voltage = 300.0u"V",
			anode_mass_flow_rate = 5.0u"mg/s",
			propellant = het.Xenon,
		)

		@test config.domain == (0.0, 0.08)
		@test config.anode_mass_flow_rate ≈  5e-6
		@test config.discharge_voltage ≈  300.0

		solution = het.run_simulation(config, simparams)
		avg = het.time_average(solution, 0.5u"ms")  # units are supported too, if Unitful or DynamicQuantities loaded
		@test length(avg.frames) == 1
		@test avg.t[end] == simparams.duration

		tenth_frame = solution[10]      # extract frame number 10 as a new Solution object
		@test length(tenth_frame.frames) == 1
		@test tenth_frame.frames[] == solution.frames[10]

		middle_800 = solution[101:900]  # extract frames 100:900
		@test length(middle_800.frames) == 800
		@test middle_800.t[1] == solution.t[101]
		@test middle_800.t[end] == solution.t[900]
	end
end

