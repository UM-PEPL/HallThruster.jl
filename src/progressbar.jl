using Printf

# define a new column type
mutable struct SimTimeColumn <: AbstractColumn
    segments::Vector
    measure::Measure
    style::String
    t::Float64
end

# constructor
SimTimeColumn(dt, style::String = "") = SimTimeColumn([], Measure(10, 1), style, dt)

Term.progress.update(col::SimTimeColumn, args...)::String = apply_style("[$(col.style)]$(@sprintf("%.2f μs", col.t * 1e6))[/$(col.style)]")

function make_progress_bar(niters, dt, config)
    if config.progress_interval > 0

        description = apply_style("[bold]Simulating...[/bold]")

        cols = [
            DescriptionColumn(description),
            Term.progress.SeparatorColumn(),
            BarColumn(),
            Term.progress.SeparatorColumn(),
            SimTimeColumn(dt, "bold")
        ]

        progress_bar = ProgressBar(;N = niters ÷ config.progress_interval, columns = cols)
    else
        progress_bar = nothing
    end
end

function stop_progress_bar!(progress_bar, params)
    if params.config.progress_interval > 0
        stop(progress_bar)
    end
end

function update_progress_bar!(progress_bar, params, t)
    progress_interval = params.config.progress_interval
    iteration = params.iteration[1]
    if progress_interval > 0 && iteration % progress_interval == 0
        for col in params.progress_bar.columns
            if col isa SimTimeColumn
                col.t = t
            end
        end
        update(progress_bar)
    end
end