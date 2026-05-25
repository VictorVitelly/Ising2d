program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  use head
  implicit none

  call cpu_time(starting)
  
  !Write thermalization history in a file and computes autocorrelation
  !call thermalize(2.5_dp)

  !Measure energy, magnetization, susceptibility, heat capacity and binder cumulant in
  !an interval of temperatures, (initial temp., final temp, n. of points between them)
  !call vary_temp(2.19_dp,2.43_dp,41)
  !call vary_temp(2.0_dp,3.0_dp,21)

  !Measure correlation function in an interval of temperatures
  !(initial temp., final temp, n. of points between them)
  call correlate(Tc,Tc,10)
  
  !call correlate2(Tc,3._dp,10)
  
  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

end program main
