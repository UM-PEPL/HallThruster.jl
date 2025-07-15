#include <cstdio>

#ifdef __cplusplus
extern "C" {
#endif

void two_zone_bohm(double *nu_anom, double *z, double *B, double *params, double channel_length, int num_cells) {
    double q_e = 1.60217664e-19;
    double m_e = 9.1093837e-31;

    for (int i = 0; i < num_cells; i++) {
        double cyclotron_freq = q_e * B[i] / m_e;
        if (z[i] < channel_length) {
            nu_anom[i] = params[0] * cyclotron_freq;
        } else {
            nu_anom[i] = params[1] * cyclotron_freq;
        }
    }
}

#ifdef __cplusplus
}
#endif