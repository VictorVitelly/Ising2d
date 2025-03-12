module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine initialize1(E,M,susc1,heat1)
    real(dp), dimension(:), intent(inout) :: E,M,susc1,heat1
      E=0._dp
      M=0._dp
      susc1=0._dp
      heat1=0._dp
  end subroutine initialize1

  subroutine initialize2(corr1,corr2)
    real(dp), dimension(N,Nmsrs), intent(inout) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(inout) :: corr2
      corr1=0._dp
      corr2=0._dp
  end subroutine initialize2

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

  subroutine correlation(spin,k,corr1,corr2)
    integer(i4), dimension(N,N), intent(in) :: spin
    integer(i4), intent(in) :: k
    !real(dp) :: M
    real(dp), dimension(N,Nmsrs), intent(inout) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(inout) :: corr2
    real(dp), dimension(N) :: spinvec
    integer(i4) :: i1,i2
    !M=0._dp
    spinvec=0._dp
    do i1=1,N
      do i2=1,N
        !M=M+phi(i1,i2)
        spinvec(i1)=spinvec(i1)+real(spin(i1,i2),dp)
      end do
    end do
    do i1=1,N
      corr1(i1,k)=spinvec(i1)
      do i2=1,N
        corr2(i1,i2,k)=spinvec(i1)*spinvec(i2)
      end do
    end do
  end subroutine correlation

  subroutine correlation_function(corr1,corr2,CF,CFprom)
    real(dp), dimension(N,Nmsrs), intent(in) :: corr1
    real(dp), dimension(N,N,Nmsrs), intent(in) :: corr2
    real(dp), dimension(N,N), intent(out) :: CF,CFprom
    real(dp), dimension(N) :: corr1prom,corr1delta
    real(dp), dimension(N,N) :: corr2prom,corr2delta
    integer(i4) :: i1,i2
    corr1prom=0._dp
    corr2prom=0._dp
    corr1delta=0._dp
    corr2delta=0._dp
    call mean_vector(corr1,corr1prom,corr1delta)
    call mean_matrix(corr2,corr2prom,corr2delta)
    do i1=1,N
      do i2=1,N
        CF(i1,i2)=corr2prom(i1,i2)-corr1prom(i1)*corr1prom(i2)
        CFprom(i1,i2)=Sqrt((corr2delta(i1,i2))**2+(corr1prom(i1)*corr1delta(i2))**2 +(corr1prom(i2)*corr1delta(i1) )**2)
      end do
    end do
  end subroutine correlation_function

end module measurements
