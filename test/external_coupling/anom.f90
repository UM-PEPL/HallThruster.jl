module Anom

use iso_fortran_env
use iso_c_binding

contains

subroutine two_zone_bohm(nu_anom, z, B, params, channel_length, num_cells) bind(c)
    implicit none
    integer(int64), intent(in):: num_cells
    real(real64), intent(out):: nu_anom(num_cells)
    real(real64), intent(in):: z(num_cells), B(num_cells), params(2), channel_length
    real(real64), parameter:: q_e = 1.60217663e-19, m_e = 9.1093837e-31
    real(real64):: cyclotron_freq = 0.0
    integer(int64):: i = 0

    do i = 1, num_cells
        cyclotron_freq = q_e * B(i) / m_e
        if (z(i) .lt. channel_length) then
            nu_anom(i) = params(1) * cyclotron_freq
        else
            nu_anom(i) = params(2) * cyclotron_freq
        endif
    end do

end subroutine two_zone_bohm

end module Anom