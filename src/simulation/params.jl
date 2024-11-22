@kwdef struct SimParams{G <: GridSpec, C <: CurrentController}
    # Grid setup
    grid::G = EvenGrid(0)

    # Timestep control
    dt::Float64 = 1e-8
    min_dt::Float64 = 1e-10
    max_dt::Float64 = 1e-7
    CFL::Float64 = 0.799
    adaptive::Bool = false
    max_small_steps::Int = 100

    # PID control
    current_control::C = NoController()

    # Reporting
    verbose::Bool = true
    show_errors::Bool = true
end
