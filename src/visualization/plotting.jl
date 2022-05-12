function load_hallis_output(output_path)
    output_headers = [
        :z, :ne, :ϕ, :Te, :Ez, :Br, :nn, :ndot, :μe, :μen, :μbohm, :μwall, :μei,
    ]
    output = DataFrame(readdlm(output_path, Float64), output_headers)
    output.ωce = output.Br * 1.6e-19 / 9.1e-31
    replace!(output.nn, 0.0 => 1e12)
    return output[1:end-1, :]
end

function load_hallis_output_2D(output_path)
    output_headers = [
        :z, :ne, :ϕ, :Te, :Ez, :Br, :nn, :ndot, :μe, :μen, :μbohm, :μwall, :μei, :f_e, :f_en, :f_bohm, :f_wall, :f_ei
    ]
    output = DataFrame(readdlm(output_path, Float64), output_headers)
    output.ωce = output.Br * 1.6e-19 / 9.1e-31
    replace!(output.nn, 0.0 => 1e12)
    #output.Te .*= 2/3
    return output[1:end-1, :]
end

function load_hallis_for_input()
    hallis = load_hallis_output("landmark/Av_PLOT_HALLIS_1D_01.out")
    ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, hallis.ϕ)
    grad_ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, -hallis.Ez)
    return ϕ_hallis, grad_ϕ_hallis
end

function load_slicezr(slicezr_path)
    slicezr_headers = [
        :z, :ne, :nn, :Te, :ϕ, :Ωe, :ni_1_1, :ni_1_2, :ni_2_1, :uiz_1_1, :uPIC, :uiz_2_1
    ]

    slicezr = DataFrame(readdlm(slicezr_path, Float64), slicezr_headers)
    slicezr.z = slicezr.z./100
    slicezr.Te = slicezr.Te .* 3/2
    return slicezr
end

function load_slice(slice_path)
    slice_headers = [
        :z, :ne, :nn, :Te, :ϕ, :f_en, :f_ei, :f_wall, :f_iz, :Ωe, :ni_1_1, :ni_2_1, :ni_1_3, :uiz_1_1, :uir_1_1
    ]

    slice = DataFrame(readdlm(slice_path, Float64), slice_headers)
    slice.f_e = slice.f_en + slice.f_ei + slice.f_wall
    slice.ωce = slice.Ωe .* slice.f_e
    slice.z = slice.z/100
    slice.Te = slice.Te .* 3/2
    slice.uiz_1_1 = slice.uiz_1_1./1000
    return slice
end

function current_spectrum(sol)

    # Number of samples
    N = length(sol.t)

    # Sample rate (Hz)
    sample_rate = 1 / (sol.t[2] - sol.t[1])

    # Get solution discharge current
    current = [sol[:Id][i][] for i in eachindex(sol.t)]

    # Fourier transform current to obtain spectrum
    amplitude = fft(current) |> fftshift .|> abs
    freqs = fftfreq(N, sample_rate) |> fftshift

    return freqs, amplitude
end