!----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------
! This module converts a dense matrix to the Sparse CSR format. It is used in the solif model, because PARDISO needs
! the matrix in this format. other routines may perform better; we do not care here
! SEE INTEL MKL DOCUMENTATION FOR DETAILS ON THE SPARSE FORMAT
!----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------


module CompressMatrixCSR
  implicit none
  contains

  !##################################################################################################################################################################
  subroutine CompressMatrix(Jacobi, values, columns, rows)
    use share_vars
    implicit none
    real (kind=pr), dimension(1:2*ns+4,1:2*ns+4), intent (in) :: Jacobi
    real (kind=pr), intent (inout) :: values(:)
    integer, intent (inout) :: columns(:)
    integer, intent (inout) :: rows(:) 
    integer :: i,j,n,k,newrow
    integer, dimension(2) :: bound

    values = 0.0
    columns = 0
    rows = 0
    n=1
    k=1
    newrow=1
    do j = 1,2*ns+4 !loop over rows
      do i = 1,2*ns+4 
	if (Jacobi(i,j).ne.0.0) then
	  if (newrow == 1) then
	    rows(k)=n
	    k = k + 1
	    newrow = 0
	  endif
	  values(n) = Jacobi(i,j)
	  columns(n) = i;
	  n = n + 1
	endif
      enddo
      newrow = 1
    enddo
    rows(k) = n

  end subroutine CompressMatrix
end module CompressMatrixCSR
