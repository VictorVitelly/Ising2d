module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use functions
  implicit none

contains

  subroutine measure(spin,E,M)
    integer(i4), dimension(:,:), intent(in) :: spin
    real(dp), intent(out) :: E, M
      E=Hamilt(spin)
      M=abs(Magnet(spin))
  end subroutine measure

end module measurements
