@recipe function f(sol::Solution, frame::Int = length(sol.u))
    (;z_cell, ionization_reactions, ncharge, L_ch) = sol.params
    model = sol.params.config.ionization_model

    subplot_width = 500
    subplot_height = 500

    plot_layout = (2, 4)
    layout := plot_layout

    #delete!(plotattributes, :frame)

    width = subplot_width * plot_layout[2]
    height = subplot_height * plot_layout[1]

    size := (width, height)

    z_normalized = z_cell ./ L_ch

    label_user = get(plotattributes, :label, "")
    label_mod = if sol.retcode == :LIF_data
        ""
    else
        label_user
    end
    label --> label_mod
    linewidth --> 2
    yscale --> :identity

    # Common options
    x := z_normalized
    xlabel := "z / L"
    margin := 10mm
    framestyle := :box

    comma = isempty(label_mod) ? "" : ", "

    # Plot neutral density
    @series begin
        y := sol[:nn][frame]
        ylabel := "Density (m⁻³)"
        subplot := 1
        title := "Neutral density"
        label := label_mod * comma * "Total"
        color := 1
        ()
    end

    ne = sol[:ne][frame]
    Tev = sol[:Tev][frame]

    charge_labels = [string(sol.params.config.propellant(Z)) for Z in 1:ncharge]

    for Z in 1:ncharge
        if sol.retcode != :LIF_data
            # Plot ion density
            @series begin
                y := sol[:ni, Z][frame]
                ylabel := "Density (m⁻³)"
                subplot := 2
                label := label_mod * comma * charge_labels[Z]
                title := "Plasma density"
                color := Z
                ()
            end
        end
    end

    if ncharge > 1
        # Plot electron density
        @series begin
            y := sol[:ne][frame]
            ylabel := "Density (m⁻³)"
            subplot := 2
            label := label_mod * comma * "ne"
            color := ncharge+1
            ()
        end
    end

    for Z in 1:ncharge
        # Plot ion velocity
        @series begin
            y := sol[:ui, Z][frame] ./ 1000
            ylabel := "Ion velocity (km/s)"
            label := label_user * comma * charge_labels[Z]
            yscale := :identity
            legend:= :bottomright
            title := "Ion velocity"
            subplot := 3
            color := Z
            ()
        end
    end

    for Z in 1:ncharge
        # Plot rate of production of all three species
        ionization_rate = fill(eps(Float64), length(z_normalized))
        if isempty(ionization_reactions)
            ionization_rate .*= NaN
        else
            for rxn in ionization_reactions
                if rxn.product.Z == Z
                    if rxn.reactant.Z == 0
                        reactant_density = sol[:nn][frame]
                    else
                        reactant_density = sol[:ni, rxn.reactant.Z][frame]
                    end
                    ionization_rate .+= reactant_density .* ne .* [rate_coeff(model, rxn, 3/2 .* T) for T in Tev]
                end
            end
        end
        @series begin
            y := ionization_rate
            ylabel := "Ionization rate (m⁻³/s)"
            subplot := 5
            label := label_mod * comma * charge_labels[Z]
            title := "Ionization rate"
            color := Z
            ()
        end
    end

    # Plot potential
    @series begin
        y := sol[:ϕ][frame]
        ylabel := "Potential (V)"
        yscale := :identity
        title := "Potential"
        subplot := 4
        label := label_mod
        color := 1
        ()
    end

    # Plot electron temperature
    @series begin
        y := sol[:Tev][frame]
        ylabel := "Electron temperature (eV)"
        yscale := :identity
        subplot := 6
        title := "Electron temperature"
        label := label_mod
        color := 1
        ()
    end

    # Plot electron velocity
    @series begin
        y := sol[:ue][frame] ./ 1000
        ylabel := "Electron velocity (km/s)"
        yscale := :identity
        subplot := 7
        title := "Cross-field electron velocity"
        label := label_mod
        color := 1
        ()
    end

    # Plot electric field
    @series begin
        y := -sol[:∇ϕ][frame]
        ylabel := "Electric field (V/m)"
        yscale := :identity
        subplot := 8
        title := "Electric field"
        label := label_mod
        color := 1
        ()
    end
end

@userplot CollisionPlot

@recipe function f(c::CollisionPlot)
    cond1 = length(c.args) > 2
    cond2 = !(typeof(c.args[1]) <: Solution)
    cond3 = (length(c.args) == 2 && !(typeof(c.args[2]) <: Integer))
    if  cond1 || cond2 || cond3
        error("Collision plot takes a Hall thruster Solution as a first argument and an optional integer as a second argument.\nGot: $(typeof(c.args))")
    end

    sol = c.args[1]

    if length(c.args) == 2
        frame = c.args[2]
    else
        frame = length(sol.u)
    end

    yaxis := :log
    xlabel := "z/L"
    ylabel := "Frequency (Hz)"
    legend := :outertop
    framestyle := :box
    ylims := get(plotattributes, :ylims, (1e5, 1e10))
    freqs = get(plotattributes, :freqs, [:ωce, :νan, :νen, :νei, :νiz, :νex, :νw])

    zs = sol.params.z_cell ./ sol.params.L_ch

    if :ωce in freqs
        @series begin
            label := "Electron cyclotron frequency"
            lw := 2
            ls := :dash
            zs, sol[:ωce][1]
        end
    end

    if :νan in freqs
        @series begin
            label := "Anomalous"
            ys = sol[:νan][frame]
            inds = findall(==(0), ys)
            ys[inds] .= NaN
            zs, ys
        end
    end

    if :νen in freqs
        @series begin
            label := "Electron-neutral"
            ys = sol[:νen][frame]
            inds = findall(==(0), ys)
            ys[inds] .= NaN
            zs, ys
        end
    end

    if :νei in freqs
        @series begin
            label := "Electron-ion"
            ys = sol[:νei][frame]
            inds = findall(==(0), ys)
            ys[inds] .= NaN
            zs, ys
        end
    end

    if :νiz in freqs
        @series begin
            label := "Ionization"
            ys = sol[:νiz][frame]
            inds = findall(==(0), ys)
            ys[inds] .= NaN
            zs, ys
        end
    end

    if :νex in freqs
        @series begin
            label := "Excitation"
            ys = sol[:νex][frame]
            inds = findall(==(0), ys)
            ys[inds] .= NaN
            zs, ys
        end
    end

    if :νw in freqs
        @series begin
            label := "Wall"
            ys = sol[:νew_momentum][frame]
            inds = findall(==(0), ys)
            ys[inds] .= NaN
            zs, ys
        end
    end

    @series begin
        lc := :black
        lw := 2
        label := "Total collision freq"
        zs, sol[:νe][frame]
    end
end
