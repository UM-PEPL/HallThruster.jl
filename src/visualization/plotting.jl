

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
