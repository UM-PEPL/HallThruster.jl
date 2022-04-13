@recipe function f(sol::HallThrusterSolution)
    (;z_cell, config, ionization_reactions) = sol.params
    (;ncharge, thruster) = config

    frame = length(sol.u)

    subplot_width = 500
    subplot_height = 500

    plot_layout = (2, 4)
    layout := plot_layout

    width = subplot_width * plot_layout[2]
    height = subplot_height * plot_layout[1]

    size := (width, height)

    z_normalized = z_cell ./ thruster.geometry.channel_length

    yscale_user = get(plotattributes, :yscale, :identity)
    label_user = get(plotattributes, :label, "")
    linewidth := get(plotattributes, :linewidth, 2)

    # Common options
    x := z_normalized
    xlabel := "z / L"
    margin := 10Measures.mm
    framestyle := :box

    # Plot neutral density
    @series begin
        y := sol[:nn][frame]
        yscale := yscale_user
        ylabel := "Density (m⁻³)"
        subplot := 1
        title := "Neutral density"
        ()
    end

    ne = sol[:ne][frame]
    Tev = sol[:Tev][frame]

    for Z in 1:ncharge
        # Plot ion density
        @series begin
            y := sol[:ni, Z][frame]
            yscale := yscale_user
            ylabel := "Density (m⁻³)"
            subplot := 2
            label := ifelse(!isempty(label_user),  label_user * ", Z = $Z", "Z = $Z")
            ()
        end

        # Plot ion velocity
        @series begin
            y := sol[:ui, Z][frame] ./ 1000
            ylabel := "Ion velocity (km/s)"
            legend := :bottomright
            label := ifelse(!isempty(label_user),  label_user * ", Z = $Z", "Z = $Z")
            yscale := :identity
            title := "Ion velocity"
            subplot := 3
            ()
        end

        # Plot rate of production of all three species
        ionization_rate = zeros(length(z_normalized))
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
        @series begin
            y := ionization_rate
            yscale := yscale_user
            ylabel := "Ionization rate (m⁻³/s)"
            subplot := 5
            label := ifelse(!isempty(label_user),  label_user * ", Z = $Z", "Z = $Z")
            title := "Ionization rate"
            ()
        end
    end

    if config.ncharge > 1
        # Plot electron density
        @series begin
            y := sol[:ne][frame]
            yscale := yscale_user
            ylabel := "Density (m⁻³)"
            subplot := 2
            label := ifelse(!isempty(label_user),  label_user * ", ne", "ne")
            title := "Plasma densities"
            ()
        end
    end

    # Plot potential
    @series begin
        y := sol[:ϕ_cell][frame]
        ylabel := "Potential (V)"
        label := label_user
        yscale := :identity
        title := "Potential"
        subplot := 4
        ()
    end

    # Plot electron temperature
    @series begin
        y := sol[:Tev][frame]
        ylabel := "Electron temperature (eV)"
        label := label_user
        yscale := :identity
        subplot := 6
        title := "Electron temperature"
        ()
    end

    # Plot electron velocity
    @series begin
        y := sol[:ue][frame] ./ 1000
        ylabel := "Electron velocity (km/s)"
        yscale := :identity
        subplot := 7
        title := "Cross-field electron velocity"
        ()
    end

    # Plot electric field
    @series begin
        y := -sol[:∇ϕ][frame]
        ylabel := "Electric field (V/m)"
        yscale := :identity
        subplot := 8
        title := "Electric field"
        ()
    end
end