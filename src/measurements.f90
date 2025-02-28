module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  implicit none

contains

  subroutine initialize1(E,M,susc1,heat1)
    real(dp), dimension(:), intent(inout) :: E,M,susc1,heat1
      E=0._dp
      M=0._dp
      susc1=0._dp
      heat1=0._dp
  end subroutine initialize1

  subroutine measure(spin,E,M,susc1,heat1)
    integer(i4), dimension(:,:), intent(in) :: spin
    real(dp), intent(out) :: E, M
    real(dp), intent(inout) :: susc1,heat1
    real(dp) :: MM
      MM=Magnet(spin)
      E=Hamilt(spin)
      M=abs(MM)
      susc1=susc1+(MM**2)
      heat1=heat1+(E**2)
  end subroutine measure

  subroutine divideN(susc1,heat1)
    real(dp), intent(inout) :: susc1,heat1
    susc1=susc1/real(Nmsrs,dp)
    heat1=heat1/real(Nmsrs,dp)
  end subroutine divideN

end module measurements
