@recipe function f(sol::Solution, frame::Int = length(sol.u))
    (;z_cell, ionization_reactions, ncharge, L_ch) = sol.params

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

    # Plot neutral density
    @series begin
        y := sol[:nn][frame]
        ylabel := "Density (m⁻³)"
        subplot := 1
        title := "Neutral density"
        label := label_mod
        ()
    end

    ne = sol[:ne][frame]
    Tev = sol[:Tev][frame]

    for Z in 1:ncharge
        if sol.retcode != :LIF_data
            # Plot ion density
            @series begin
                y := sol[:ni, Z][frame]
                ylabel := "Density (m⁻³)"
                subplot := 2
                label := ifelse(!isempty(label_mod),  label_mod * ", Z = $Z", "")
                title := "Plasma density"
                ()
            end
        end

        # Plot ion velocity
        @series begin
            y := sol[:ui, Z][frame] ./ 1000
            ylabel := "Ion velocity (km/s)"
            legend := :bottomright
            label := ifelse(!isempty(label_user),  label_user * ", Z = $Z", "")
            yscale := :identity
            title := "Ion velocity"
            subplot := 3
            ()
        end

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
                    ionization_rate .+= reactant_density .* ne .* rxn.rate_coeff.(3/2 .* Tev)
                end
            end
        end
        @series begin
            y := ionization_rate
            ylabel := "Ionization rate (m⁻³/s)"
            subplot := 5
            label := ifelse(!isempty(label_mod), label_mod * ", Z = $Z", "")
            title := "Ionization rate"
            ()
        end
    end

    if ncharge > 1
        # Plot electron density
        @series begin
            y := sol[:ne][frame]
            ylabel := "Density (m⁻³)"
            subplot := 2
            label := ifelse(!isempty(label_mod),  label_mod * ", ne", "ne")
            ()
        end
    end

    # Plot potential
    @series begin
        y := sol[:ϕ_cell][frame]
        ylabel := "Potential (V)"
        yscale := :identity
        title := "Potential"
        subplot := 4
        label := label_mod
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
        ()
    end
end