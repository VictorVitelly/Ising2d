module statistics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  implicit none

contains

  subroutine cold_start(spin)
    integer(i4), dimension(:,:), intent(out) :: spin
    spin=-1
  end subroutine cold_start

  subroutine hot_start(spin)
    integer(i4), dimension(:,:), intent(out) :: spin
    integer(i4) :: i,j
    real(dp) :: r1
    do i=1,size(spin,dim=1)
      do j=1,size(spin,dim=2)
        call random_number(r1)
        if(r1 .le. 0.5_dp ) then
          spin(i,j)=1
        else
          spin(i,j)=-1
        end if
      end do
    end do
  end subroutine hot_start

  subroutine montecarlo(spin,T)
    integer(i4), dimension(:,:), intent(inout) :: spin
    real(dp), intent(in) :: T
    integer(i4) :: i1,i2
    real(dp) :: dH,r1,p
    do i1=1,size(spin,dim=1)
      do i2=1,size(spin,dim=2)
        dH=DeltaH(spin,i1,i2)
        if(dH .le. 0._dp) then
          spin(i1,i2)=-1*spin(i1,i2)
        else
          call random_number(r1)
          p=Exp(-dH/T)
          if(r1 < p ) then
            spin(i1,i2)=-1*spin(i1,i2)
          end if
        end if
      end do
    end do
  end subroutine montecarlo

  !Error statistics

  subroutine standard_error(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: variance
    integer(i4) :: k,N
    N=size(x)
    deltay=0._dp
    variance=0._dp
    do k=1,N
      variance=variance+(x(k) -y)**2
    end do
    variance=variance/real(N-1,dp)
    deltay=Sqrt(variance/real(N,dp))
  end subroutine standard_error

  subroutine jackknife(x,Mbins,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    integer(i4), intent(in) :: Mbins
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: jackk
    real(dp),dimension(Mbins) :: xmean
    integer(i4) :: k,N,i
      N=size(x)
      deltay=0._dp
      jackk=0._dp
      xmean=0._dp
      do i=1,Mbins
        do k=1,N
          if(k .le. (i-1)*N/Mbins) then
            xmean(i)=xmean(i)+x(k)
          else if(k > i*N/Mbins) then
            xmean(i)=xmean(i)+x(k)
          end if
        end do
        xmean(i)=xmean(i)/(real(N,dp) -real(N/Mbins,dp))
      end do
      do k=1,Mbins
        jackk=jackk+(xmean(k)-y )**2
      end do
      deltay=Sqrt(real(Mbins-1,dp)*jackk/real(Mbins,dp))
  end subroutine jackknife

  subroutine mean_scalar(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y,deltay
    integer(i4) :: k,N
    N=size(x)
    y=0._dp
    do k=1,N
      y=y+x(k)
    end do
    y=y/real(N,dp)
    call standard_error(x,y,deltay)
    !call jackknife(x,10,y,deltay)
  end subroutine mean_scalar

end module statistics
