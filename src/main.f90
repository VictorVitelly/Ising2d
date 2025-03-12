program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  !call thermalize(2.25_dp)

  open(20, file = 'data/energy.dat', status = 'replace')
  open(30, file = 'data/magnetization.dat', status = 'replace')
  open(40, file = 'data/susceptibility.dat', status = 'replace')
  open(50, file = 'data/heat.dat', status = 'replace')
  allocate(spin(N,N))
  allocate(E(Nmsrs))
  allocate(M(Nmsrs))
  allocate(susc1(Nmsrs))
  allocate(heat1(Nmsrs))

  call correlate(3._dp)

  do j=0,30
    write(*,*) j
    call initialize1(E,M,susc1,heat1)
    call cold_start(spin)
    k=0
    T=1._dp+0.1_dp*real(j,dp)
    !T=2.24_dp+0.005_dp*real(j,dp)
    do i=1,sweeps
      call montecarlo(spin,T )
      if(i>thermalization .and. mod(i,eachsweep)==0) then
        k=k+1
        call measure(spin,E(k),M(k),susc1(k),heat1(k))
      end if
    end do
    call mean_scalar(E,Emean,deltaE)
    call mean_scalar(M,Mmean,deltaM)
    call heat_jackk(susc1,M,suscmean,deltasusc)
    call heat_jackk(heat1,E,heatmean,deltaheat)
    write(20,*) T, Emean/(real(N**2,dp) ), deltaE/(real(N**2,dp) )
    write(30,*) T, Mmean/(real(N**2,dp) ), deltaM/(real(N**2,dp) )
    write(40,*) T, (suscmean)/(real(N**2,dp)*T ), deltasusc/(real(N**2,dp)*T )
    write(50,*) T, (heatmean)/(real(N**2,dp)*(T**2) ), deltaheat/(real(N**2,dp)*(T**2) )
  end do

  close(20)
  close(30)
  close(40)
  close(50)
  deallocate(spin,E,M)
  deallocate(susc1,heat1)

contains

  subroutine thermalize(T)
  real(dp), intent(in) :: T
  integer(i4) :: i
  integer(i4), allocatable :: spin(:,:)
  open(10, file = 'data/therm.dat', status = 'replace')
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

  subroutine correlate(T)
  real(dp), intent(in) :: T
  integer(i4) :: i,k
  integer(i4), allocatable :: spin(:,:)
  real(dp), allocatable :: corr1(:,:)
  real(dp), allocatable :: corr2(:,:,:)
  real(dp), allocatable :: CF(:,:),CFprom(:,:)
  open(60, file = 'data/corrfunc.dat', status = 'replace')
    allocate(corr1(N,Nmsrs))
    allocate(corr2(N,N,Nmsrs))
    allocate(CF(N,N))
    allocate(CFprom(N,N))
    allocate(spin(N,N))
    call cold_start(spin)
    call initialize2(corr1,corr2)
    do i=1,sweeps
      call montecarlo(spin,T)
      if(i>thermalization .and. mod(i,eachsweep)==0) then
        k=k+1
        call correlation(spin,k,corr1,corr2)
      end if
    end do
    call correlation_function(corr1,corr2,CF,CFprom)
    do i=1,N+1
      write(60,*) abs(i-1), CF(iv(i),1), CFprom(iv(i),1)
    end do
    close(60)
    deallocate(spin)
    deallocate(corr1,corr2,CF,CFprom)
  end subroutine correlate


end program main
