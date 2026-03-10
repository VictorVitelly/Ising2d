module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: N=32,thermalization=2000,eachsweep=10,Nmsrs=200
    integer(i4),parameter :: Nauto=15000,Nmsrs2=120,Mbin(4)=(/5,10,15,20/)
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real :: starting,ending
    real(dp), parameter :: PI=4._dp*Atan(1.0_dp)
  
end module parameters
