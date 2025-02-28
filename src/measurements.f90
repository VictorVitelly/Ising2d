module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use functions
  implicit none

contains

  subroutine measure(spin,E,M,susc1,susc2,heat1,heat2)
    integer(i4), dimension(:,:), intent(in) :: spin
    real(dp), intent(out) :: E, M
    real(dp), intent(inout) :: susc1,susc2,heat1,heat2
    real(dp) :: MM
      MM=Magnet(spin)
      E=Hamilt(spin)
      M=abs(MM)
      susc1=susc1+(MM**2)
      susc2=susc2+M
  end subroutine measure

end module measurements
