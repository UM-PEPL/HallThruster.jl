import numpy as np
import scipy.stats.qmc as qmc

q_e = 1.6e-19
m_e = 9.1e-31

def compute_rate_coefficients(temperatures_eV, energy_eV, sigma_m2, num_samples=1024):
    """Integrate cross sections over maxwellian VDF to obtain reaction rate coefficients
    Given electron energy in eV (`energy_eV`) and reaction cross sections at those energies (`sigma_m2`),
    this function computes the reaction rate coefficient $k(T_e)$ for maxwellian electrons at
    a provided list of electron temperatures `temperatures_eV`.
    The rate coefficient is given by
    $$
    k(T_e) = \int \sigma(E) E dE
    $$
    where the energy $E$ is drawn from a maxwellian distribution function with zero speed and temperature $T_e$.
    We solve this using a quasi-monte carlo approach, by drawing a large number of low-discrepancy samples from
    the appropriate distribution and obtaining the average of $\sigma(E) E$.
    """
    thermal_speed_scale = np.sqrt(q_e / m_e)
    k = np.zeros(temperatures_eV.size)

    # obtain low-discrepancy samples of normal dist
    dist = qmc.MultivariateNormalQMC(np.zeros(3), np.eye(3))
    v = dist.random(num_samples)

    for i, T in enumerate(temperatures_eV):
        # scale velocities to proper temperature
        # compute energies corresponding to each sampled velocity vector
        speed_squared = (v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2) * T
        e = 0.5 * speed_squared
        speed = np.sqrt(speed_squared) * thermal_speed_scale
        # get cross section by interpolating on table
        sigma = np.interp(e, energy_eV, sigma_m2, left=0)
        k[i] = np.mean(sigma * speed)
    return k

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Rate coefficient calculator")

    parser.add_argument(
        "cross_section_file",
        type=str,
        help="""File containing a delimited table of electron energies and corresponding cross-sections.
        The delimiter can be set using the -d or --delimiter flags and defaults to a comma."""
    )

    parser.add_argument(
        "--output-file", "-o",
        type=str,
        help="File to write output rate coefficients to",
        default = "_rate_coeffs.csv",
    )

    parser.add_argument(
        "--delimiter", "-d",
        type=str,
        help="""Delimiter for reading and writing rate coefficient tables. Defaults to a comma.
        If seperate input and output delimiters are desired, the --input-delimiter and --output-delimiter flags are available""",
        default=",",
    )

    parser.add_argument(
        "--input-delimiter",
        type=str,
        help="""Expected delimiter for the input cross-section file. Defaults to the value of --delimiter if not set.""",
    )

    parser.add_argument(
        "--output-delimiter",
        type=str,
        help="""Expected delimiter for the output rate coefficient file. Defaults to the value of --delimiter if not set.""",
    )

    args = parser.parse_args()

    if args.input_delimiter is not None:
        input_delimiter = args.input_delimiter
    else:
        input_delimiter = args.delimiter

    if args.output_delimiter is not None:
        output_delimiter = args.output_delimiter
    else:
        output_delimiter = args.delimiter

    if input_delimiter == "\\t":
        input_delimiter = "\t"
    if output_delimiter == "\\t":
        output_delimiter = "\t"

    table = np.genfromtxt(args.cross_section_file, dtype=np.float64, delimiter=input_delimiter)
    energies = np.concatenate((np.arange(0,10,0.5), np.arange(10, 30, 2), np.arange(30, 60, 5), np.arange(60, 160, 10)))
    temps = energies/1.5
    k = compute_rate_coefficients(temps, table[:, 0], table[:, 1])

    header = "Mean energy [eV]" + output_delimiter + "Rate coefficient [m^3/s]"
    rate_coeff_table = np.stack((energies, k), axis=-1)
    np.savetxt(args.output_file, rate_coeff_table, header=header,comments="",delimiter=output_delimiter, fmt=["%.1f", "%.12e"])
