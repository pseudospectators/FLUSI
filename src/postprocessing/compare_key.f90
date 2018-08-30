!-------------------------------------------------------------------------------
! ./flusi --postprocess --compare-keys mask_00000.key saved.key
!-------------------------------------------------------------------------------
! compares to *.key files if they're equal
subroutine compare_key(key1,key2)
  use vars
  use mpi
  use module_helpers, only : check_file_exists
  implicit none
  character(len=*), intent(in) :: key1,key2
  real(kind=pr) :: a1,a2,b1,b2,c1,c2,d1,d2,t1,t2,q1,q2
  real(kind=pr) :: e1,e2,e3,e4,e0,e5
  integer ::mpicode

  call check_file_exists(key1)
  call check_file_exists(key2)

  open  (14, file = key1, status = 'unknown', action='read')
  read (14,'(6(es15.8,1x))') t1,a1,b1,c1,d1,q1
  close (14)

  open  (14, file = key2, status = 'unknown', action='read')
  read (14,'(6(es15.8,1x))') t2,a2,b2,c2,d2,q2
  close (14)

  write (*,'("present  : time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
  t1,a1,b1,c1,d1,q1

  write (*,'("reference: time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
  t2,a2,b2,c2,d2,q2

  ! errors:
  if (dabs(t2)>=1.0d-7) then
    e0 = dabs( (t2-t1) / t2 )
  else
    e0 = dabs( (t2-t1) )
  endif

  if (dabs(a2)>=1.0d-7) then
    e1 = dabs( (a2-a1) / a2 )
  else
    e1 = dabs( (a2-a1) )
  endif

  if (dabs(b2)>=1.0d-7) then
    e2 = dabs( (b2-b1) / b2 )
  else
    e2 = dabs( (b2-b1) )
  endif

  if (dabs(c2)>=1.0d-7) then
    e3 = dabs( (c2-c1) / c2 )
  else
    e3 = dabs( (c2-c1) )
  endif

  if (dabs(d2)>=1.0d-7) then
    e4 = dabs( (d2-d1) / d2 )
  else
    e4 = dabs( (d2-d1) )
  endif

  if (dabs(q2)>=1.0d-7) then
    e5 = dabs( (q2-q1) / q2 )
  else
    e5 = dabs( (q2-q1) )
  endif

  write (*,'("err(rel) : time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
  e0,e1,e2,e3,e4,e5

  if ((e1<1.d-4) .and. (e2<1.d-4) .and. (e3<1.d-4) .and. (e4<1.d-4) .and. (e0<1.d-4) .and. (e5<1.d-4)) then
    ! all cool
    write (*,*) "OKAY..."

    ! on some machines, returning an exit code (exit(1)) does not work
    ! so write your exit code in a small txt file as well. this allows unit tests
    ! on turing.
    if (root) then
      open (15, file='return', status='replace')
      write(15,'(i1)') 0
      close(15)
    endif
    call MPI_FINALIZE(mpicode)
    call exit(0)
  else
    ! very bad
    write (*,*) "ERROR"

    ! on some machines, returning an exit code (exit(1)) does not work
    ! so write your exit code in a small txt file as well. this allows unit tests
    ! on turing.
    if (root) then
      open (15, file='return', status='replace')
      write(15,'(i1)') 1
      close(15)
    endif
    call MPI_FINALIZE(mpicode)
    call exit(1)
  endif
end subroutine compare_key
