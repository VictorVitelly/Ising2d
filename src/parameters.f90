module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: N=16,thermalization=1000,eachsweep=100,Nmsrs=10000
    integer(i4),parameter :: Mbins=10
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real :: starting,ending

end module parameters
