program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  !call thermalize(1.5_dp)

  open(20, file = 'data/energy.dat', status = 'replace')
  open(30, file = 'data/magnetization.dat', status = 'replace')
  allocate(spin(N,N))
  allocate(E(Nmsrs))
  allocate(M(Nmsrs))

  do j=0,30
    call cold_start(spin)
    k=0
    T=1._dp+0.1_dp*real(j,dp)
    do i=1,sweeps
      call montecarlo(spin,T )
      if(i>thermalization .and. mod(i,eachsweep)==0) then
        k=k+1
        call measure(spin,E(k),M(k))
      end if
    end do
    call mean_scalar(E,Emean,deltaE)
    call mean_scalar(M,Mmean,deltaM)
    write(20,*) T, Emean/(real(N**2,dp) ), deltaE/(real(N**2,dp) )
    write(30,*) T, Mmean/(real(N**2,dp) ), deltaM/(real(N**2,dp) )
  end do

  close(20)
  close(30)
  deallocate(spin,E,M)

contains

  subroutine thermalize(T)
  real(dp), intent(in) :: T
  integer(i4) :: i
  integer(i4), allocatable :: spin(:,:)
  open(10, file = 'data/therm_t1p5.dat', status = 'replace')
  allocate(spin(N,N))
    call cold_start(spin)
    do i=1,1000
      if(i==1 .or. mod(i,2)==0 ) then
        write(10,*) i, Hamilt(spin)/(real(N**2,dp) )
      end if
      call montecarlo(spin,T )
    end do
  close(10)
  deallocate(spin)
  end subroutine thermalize


end program main
