module functions
    use iso_fortran_env, only : dp => real64, i4 => int32
    use parameters, only : N
    implicit none

contains

  function iv(i)
    integer(i4), intent(in) :: i
    integer(i4) :: iv
    if(i==N+1) then
      iv=1
    else if(i==0) then
      iv=N
    else
      iv=i
    end if
  end function iv

  function Hamilt(spin)
    integer(i4), dimension(:,:), intent(in) :: spin
    real(dp) :: Hamilt,neigh
    integer(i4) :: i,j
    Hamilt=0._dp
    do i=1,size(spin,dim=1)
      do j=1,size(spin,dim=2)
        !neigh=real(spin(iv(i+1),j)+spin(iv(i-1),j)+spin(i,iv(j+1))+spin(i,iv(j-1)),dp)
        !Hamilt=Hamilt-real(spin(i,j),dp)*neigh/2._dp
        neigh=real(spin(iv(i+1),j)+spin(i,iv(j+1)),dp)
        Hamilt=Hamilt-real(spin(i,j),dp)*neigh
      end do
    end do
  end function Hamilt

  function DeltaH(spin,i,j)
    integer(i4), dimension(:,:), intent(in) :: spin
    integer(i4),intent(in) :: i,j
    real(dp) :: DeltaH,neigh
    neigh=real(spin(iv(i+1),j)+spin(iv(i-1),j)+spin(i,iv(j+1))+spin(i,iv(j-1)),dp)
    DeltaH=2._dp*real(spin(i,j),dp)*neigh
  end function DeltaH


end module functions
