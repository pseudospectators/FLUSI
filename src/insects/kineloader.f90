subroutine load_kine_init(mpirank)
  use mpi
  use kine
  implicit none

  integer, intent(in) :: mpirank
  integer :: j, mpicode
  real (kind=prk) :: t_period, r_wing 

  if (mpirank==0) then
    open (10, file = 'data_kin.dat', form='formatted', status='old') 
    read (10, *) t_period ! stroke period in s, for normalization
    read (10, *) r_wing   ! wing length in mm, for normalization
    read (10, *) nk
  endif

  call MPI_BCAST(nk,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode)
  allocate(vec_t(nk))
  allocate(vec_phi(nk))
  allocate(vec_phi_dt(nk))
  allocate(vec_alpha(nk))
  allocate(vec_alpha_dt(nk))
  allocate(vec_theta(nk))
  allocate(vec_theta_dt(nk))
  allocate(vec_pitch(nk))
  allocate(vec_pitch_dt(nk))
  allocate(vec_vert(nk))
  allocate(vec_vert_dt(nk))
  allocate(vec_horz(nk))
  allocate(vec_horz_dt(nk))
  
  if (mpirank==0) then 
    ! read from file
    do j = 1, nk
      read (10, *) vec_t(j), &
        vec_phi(j),vec_alpha(j),vec_theta(j),vec_pitch(j),vec_vert(j),vec_horz(j),  &
        vec_phi_dt(j),vec_alpha_dt(j),vec_theta_dt(j),vec_pitch_dt(j),vec_vert_dt(j),vec_horz_dt(j)
    enddo
    close (10)
    print *, "load_kine_init: data read from file, nk=", nk
    ! non-dimensionalize
    vec_t(:) = vec_t(:) / t_period
    vec_vert(:) = vec_vert(:) / r_wing
    vec_horz(:) = vec_horz(:) / r_wing
    vec_phi_dt(:) = vec_phi_dt(:) * t_period
    vec_alpha_dt(:) = vec_alpha_dt(:) * t_period
    vec_theta_dt(:) = vec_theta_dt(:) * t_period
    vec_pitch_dt(:) = vec_pitch_dt(:) * t_period
    vec_vert_dt(:) = vec_vert_dt(:) * t_period / r_wing
    vec_horz_dt(:) = vec_horz_dt(:) * t_period / r_wing 
  endif
  
  call MPI_BCAST( vec_t, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_phi, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_alpha, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_theta, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_pitch, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_vert, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_horz, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_phi_dt, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_alpha_dt, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_theta_dt, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_pitch_dt, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_vert_dt, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( vec_horz_dt, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
end subroutine


subroutine load_kine_clean
  use kine
  implicit none

  deallocate(vec_t)
  deallocate(vec_phi)
  deallocate(vec_phi_dt)
  deallocate(vec_alpha)
  deallocate(vec_alpha_dt)
  deallocate(vec_theta)
  deallocate(vec_theta_dt)
  deallocate(vec_pitch)
  deallocate(vec_pitch_dt)
  deallocate(vec_vert)
  deallocate(vec_vert_dt)
  deallocate(vec_horz)

end subroutine


subroutine wing_kine_interp(t_i,phi_i,alpha_i,theta_i,phi_dt_i,alpha_dt_i,theta_dt_i)
  use kine
  implicit none

  real(kind=prk), intent(in) :: t_i
  real(kind=prk), intent(out) :: phi_i,alpha_i,theta_i,phi_dt_i,alpha_dt_i,theta_dt_i

  call hermite1d(nk,vec_t,vec_phi,vec_phi_dt,t_i,phi_i,phi_dt_i)
  call hermite1d(nk,vec_t,vec_alpha,vec_alpha_dt,t_i,alpha_i,alpha_dt_i)
  call hermite1d(nk,vec_t,vec_theta,vec_theta_dt,t_i,theta_i,theta_dt_i)
end subroutine


subroutine body_kine_interp(t_i,pitch_i,vert_i,horz_i,pitch_dt_i,vert_dt_i,horz_dt_i)
  use kine
  implicit none

  real(kind=prk), intent(in) :: t_i
  real(kind=prk), intent(out) :: pitch_i,vert_i,horz_i,pitch_dt_i,vert_dt_i,horz_dt_i

  call hermite1d(nk,vec_t,vec_pitch,vec_pitch_dt,t_i,pitch_i,pitch_dt_i)
  call hermite1d(nk,vec_t,vec_vert,vec_vert_dt,t_i,vert_i,vert_dt_i)
  call hermite1d(nk,vec_t,vec_horz,vec_horz_dt,t_i,horz_i,horz_dt_i)
end subroutine


subroutine hermite1d(n, xphi, phi, dpdx, xi, phi_interp, dpdx_interp)
  implicit none

  integer, parameter :: pr = 8
  integer :: i0,i1
  integer, intent(in) :: n
  real(kind=pr) :: x,z,ap0,ap1,ax0,ax1,dx,d2pdx2_i0,d2pdx2_i1
  real(kind=pr), intent(in) :: xi
  real(kind=pr), intent(in) :: xphi(1:n),phi(1:n),dpdx(1:n)
  real(kind=pr), intent(out) :: phi_interp,dpdx_interp

  dx = xphi(2) - xphi(1)

  i0 = floor(xi/dx)+1
  i1 = i0+1

  if ((xi<xphi(i0)).or.(xi>xphi(i1))) then
     print *, "hermite1d: not a uniform grid"
     write(*,'("xi=",es12.4," dx=",es12.4)') xi,dx
     write(*,'("xphi(i0)=",es12.4," xphi(i1)=",es12.4)') xphi(i0),xphi(i1)
     call abort()
  endif

  x = (xi-xphi(i0))/dx
  z = 1.0d0-x

  ! Phi
  ap0 = 2.0d0*x*x*x-3.0d0*x*x+1.0d0 !f(x)
  ap1 = 2.0d0*z*z*z-3.0d0*z*z+1.0d0 ! f(1.0d0-x)
  ax0 = x*x*x-2.0d0*x*x+x !g(x)
  ax1 = -(z*z*z-2.0d0*z*z+z) !-g(1.0d0-x)
  phi_interp=(phi(i0)*ap0+phi(i1)*ap1)+dx*(dpdx(i0)*ax0+dpdx(i1)*ax1)

  ! Phi_x
  if ((i0>1).and.(i1<n)) then
     d2pdx2_i0 = 0.5d0*(dpdx(i1)-dpdx(i0-1))/dx
     d2pdx2_i1 = 0.5d0*(dpdx(i1+1)-dpdx(i0))/dx
     dpdx_interp=(dpdx(i0)*ap0+dpdx(i1)*ap1)+dx*(d2pdx2_i0*ax0+d2pdx2_i1*ax1)
  else
     ap0 = 6.0d0*x*x-6.0d0*x !df(x)
     ap1 = -(6.0d0*z*z-6.0d0*z) !-df(1.0d0-x)
     ax0 = 3.0d0*x*x-4.0d0*x+1.0d0 !dg(x)
     ax1 = 3.0d0*z*z-4.0d0*z+1.0d0 ! dg(1.0d0-x)
     dpdx_interp=(phi(i0)*ap0+phi(i1)*ap1)/dx+(dpdx(i0)*ax0+dpdx(i1)*ax1)
  endif

end subroutine


