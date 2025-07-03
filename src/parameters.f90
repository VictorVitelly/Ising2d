module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: N=32,thermalization=1000,eachsweep=100,Nmsrs=10000
    integer(i4),parameter :: Mbins=10,Nauto=10000,Nmsrs2=1500,Mbin(3)=(/5,10,20/)
    real(dp), parameter :: q=1.2_dp
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real :: starting,ending

end module parameters
