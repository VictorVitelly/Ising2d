module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine initialize(corr1,corr2)
    real(dp), dimension(N), intent(inout) :: corr1
    real(dp), dimension(N,N), intent(inout) :: corr2
      corr1=0._dp
      corr2=0._dp
  end subroutine initialize

  subroutine correlation(spin,corr1,corr2)
    integer(i4), dimension(N,N), intent(in) :: spin
    real(dp), dimension(N), intent(inout) :: corr1
    real(dp), dimension(N,N), intent(inout) :: corr2
    real(dp):: spinvec(N)
    integer(i4) :: i1,i2
    spinvec=0._dp
    !do i1=1,N
    !  do i2=1,N
    !    spinvec(i1)=spinvec(i1)+real(spin(i1,i2),dp)
    !  end do
    !end do
    do i1=1,N
      do i2=1,N
        corr2(i1,i2)=corr2(i1,i2)+spin(i1,i2)*spin(1,1)
      end do
    end do
  end subroutine correlation
  
  subroutine correlationoptim(spin,corr1)
    integer(i4), dimension(N,N), intent(in) :: spin
    real(dp), dimension(N), intent(inout) :: corr1
    real(dp):: spinvec(N)
    integer(i4) :: i1,i2
    spinvec=0._dp
    do i1=1,N
      do i2=1,N
        spinvec(i1)=spinvec(i1)+real(spin(i1,i2),dp)
      end do
    end do
    do i1=1,N
        corr1(i1)=corr1(i1)+spinvec(i1)*spinvec(1)
        !corr1(i1)=corr1(i1)+spin(1,1)*spin(i1,1)
    end do
  end subroutine correlationoptim

  subroutine correlation_function(corr1,corr2,CF)
    real(dp), dimension(N), intent(in) :: corr1
    real(dp), dimension(N,N), intent(in) :: corr2
    real(dp), dimension(N), intent(out) :: CF
    integer(i4) :: i1
    do i1=1,N
      CF(i1)=corr2(i1,1)
    end do
    !CF(:)=CF(:)/real(N**2,dp)
  end subroutine correlation_function
  
  subroutine secondmomentum(CF,xi2_ave,xi2_err)
  real(dp),dimension(N,Nmsrs2),intent(in) :: CF
  real(dp),intent(out) :: xi2_ave,xi2_err 
  integer(i4) :: i1,i2
  real(dp) :: xi2(Nmsrs2),F1(Nmsrs2),F2(Nmsrs2),F12(Nmsrs2),F12_ave,F12_err
  F1(:)=0._dp
  F2(:)=0._dp
  do i1=1,Nmsrs2
    do i2=1,N
      F1(i1)=F1(i1)+CF(i2,i1)
      F2(i1)=F2(i1)+CF(i2,i1)*COS(real(i2-1,dp)*2._dp*PI/real(N,dp))
    end do
  end do
  do i1=1,Nmsrs2
    F12(i1)=F1(i1)/F2(i1)
  end do
  call mean_scalar(F12,F12_ave,F12_err)
  xi2_ave=sqrt( (F12_ave -1._dp))/(2._dp*abs(SIN(PI/real(N,dp))) ) 
  xi2_err=F12_err/(4._dp*sqrt(F12_ave-1._dp)*abs(SIN(PI/real(N,dp))) )
  write(*,*) xi2_ave,xi2_err
  end subroutine secondmomentum
  
  subroutine secondmomentum2(CF,xi2_ave,xi2_err)
  real(dp),dimension(N,N,Nmsrs2),intent(in) :: CF
  real(dp),intent(out) :: xi2_ave,xi2_err 
  integer(i4) :: i1,i2,i3
  real(dp) :: xi2(Nmsrs2),F1(Nmsrs2),F2(Nmsrs2),F12(Nmsrs2),F12_ave,F12_err
  F1(:)=0._dp
  F2(:)=0._dp
  do i1=1,Nmsrs2
    do i2=1,N
      do i3=1,N
        !F1(i1)=F1(i1)+CF(iv(i3),iv(i2),i1)
        !F2(i1)=F2(i1)+CF(iv(i3),iv(i2),i1)*COS(real(i2-1,dp)*2._dp*PI/real(N,dp))
        F1(i1)=F1(i1)+CF(i3,i2,i1)
        F2(i1)=F2(i1)+CF(i3,i2,i1)*COS(real(i2-1,dp)*2._dp*PI/real(N,dp))
      end do
    end do
  end do
  do i1=1,Nmsrs2
    F12(i1)=F1(i1)/F2(i1)
  end do
  call mean_scalar(F12,F12_ave,F12_err)
  xi2_ave=sqrt( (F12_ave -1._dp))/(2._dp*abs(SIN(PI/real(N,dp))) ) 
  xi2_err=F12_err/(4._dp*sqrt(F12_ave-1._dp)*abs(SIN(PI/real(N,dp))) )
  write(*,*) xi2_ave,xi2_err
  end subroutine secondmomentum2
  

end module measurements
