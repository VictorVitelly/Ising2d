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
    do i1=1,N
      do i2=1,N
        spinvec(i1)=spinvec(i1)+real(spin(i1,i2),dp)
      end do
    end do
    do i1=1,N
      corr1(i1)=corr1(i1)+spinvec(i1)
      do i2=1,N
        corr2(i1,i2)=corr2(i1,i2)+spinvec(i1)*spinvec(i2)
      end do
    end do
  end subroutine correlation
  
  subroutine correlationb(spin,corr1,corr2)
    integer(i4), dimension(N,N), intent(in) :: spin
    real(dp), dimension(N), intent(inout) :: corr1
    real(dp), dimension(N,N), intent(inout) :: corr2
    integer(i4) :: i1,i2
    do i1=1,N
      corr1(i1)=spin(i1,1)
      do i2=1,N
        corr2(i1,i2)=corr2(i1,i2)+spin(i1,1)*spin(i2,1)
      end do
    end do
  end subroutine correlationb

  subroutine correlation_function(corr1,corr2,CF)
    real(dp), dimension(N), intent(in) :: corr1
    real(dp), dimension(N,N), intent(in) :: corr2
    real(dp), dimension(N), intent(out) :: CF
    integer(i4) :: i1
    do i1=1,N
      CF(i1)=corr2(i1,1)-(corr1(1)**2)
    end do
    !CF(:)=CF(:)/real(N**2,dp)
  end subroutine correlation_function
  
  subroutine correlation2(CF_ave,CF_err,xi2_ave,xi2_err)
    real(dp), intent(in) :: CF_ave(N),CF_err(N)
    real(dp),intent(out) :: xi2_ave,xi2_err
    real(dp) :: F1,F2,DF1,DF2,DFTOT
    integer(i4) :: i1,j1
    xi2_ave=0._dp
    xi2_err=0._dp
    F1=0._dp
    F2=0._dp
    DF1=0._dp
    DF2=0._dp
    do j1=1,N
        F1=F1+CF_ave(j1)
        F2=F2+CF_ave(j1)*COS(real(j1-1,dp)*2._dp*PI/real(N,dp))
        DF1=DF1+(CF_ave(j1) *CF_err(j1))**2
        DF2=DF2+(CF_ave(j1)*COS(real(j1-1,dp)*2._dp*PI/real(N,dp)) *CF_err(j1) )**2
        !write(*,*) 'CF', CF_ave(j1),CF_ave(j1)*COS(real(j1,dp)*2._dp*PI/real(N,dp))
    end do
    xi2_ave=sqrt(F1/F2-1._dp)/(2._dp*SIN(PI/N))
    DF1=SQRT(DF1)
    DF2=SQRT(DF2)
    DFTOT=SQRT((DF1/F2)**2+(DF2*F1/(F2**2))**2 )
    xi2_err=DFTOT/(4._dp*sqrt(F1/F2-1._dp)*SIN(PI/N) )
    write(*,*) 'F1/F2(i)=', F1,F2, 'xi2(i)=',xi2_ave
  end subroutine correlation2
  
  subroutine autocorrelation(T,tmax,spin)
    integer(i4), intent(in) :: tmax
    real(dp), intent(in) :: T
    integer(i4), dimension(N,N), intent(inout) :: spin
    real(dp), dimension(tmax+1) :: auto,auto_delta
    real(dp) :: E(Nmsrs+tmax), auto1(Nmsrs)
    real(dp) :: E_ave,auto1_ave,autoj(tmax+1,Nauto)
    integer(i4) :: i,j,tt
    open(70, file = 'data/autocorr.dat', status = 'replace')
    do j=1,Nauto
      do i=1,Nmsrs+tmax
        !call montecarlo(spin,T)
        call cluster(spin,T)
        E(i)=Hamilt(spin)/(N**2)
      end do
      call mean_0(E,E_ave )
      
      do tt=0,tmax
        do i=1,Nmsrs
          auto1(i)=E(i)*E(i+tt)
        end do
        call mean_0(auto1,auto1_ave)
        auto=auto1_ave-(E_ave**2)
        autoj(tt+1,j)=auto1_ave-(E_ave**2)
      end do
    end do
    do tt=0,tmax
      call mean_scalar(autoj(tt+1,:),auto(tt+1),auto_delta(tt+1))
      write(70,*) tt,auto(tt+1),auto_delta(tt+1)
    end do
    close(70)
  end subroutine autocorrelation

end module measurements
