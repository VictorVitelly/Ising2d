module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: N=8,thermalization=1000,eachsweep=10,Nmsrs=10000
    integer(i4),parameter :: Mbins=10
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs

    integer(i4) :: i,k,j
    real(dp),allocatable :: E(:),M(:)
    real(dp) :: T,Emean,deltaE,Mmean,deltaM
    real(dp),allocatable :: susc1(:),heat1(:)
    real(dp) :: suscmean,deltasusc,heatmean,deltaheat


end module parameters
