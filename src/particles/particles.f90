! Variables for particles part
module particles
  ! Benjamin 03 July 2019
  ! Module to compute particles trajectories and Lagrangian statitics
  implicit none
  integer,save :: use_lagrangian
  integer,save :: reclc ! reclc is the number of bytes for I/O access direct
  integer,save :: npart ! total number of particles
  integer,save :: nper_part !number of particles allowed by proc in percentage: +nper_part% of npart/nproc
  integer,save :: npart_proc ! number of particles allowed by proc
  integer,save :: n_interp ! number of field to interpolate for lagrangian
  integer,save :: n_lagr ! number of field to save for lagrangian
  integer,save :: ilagr_acc ! compute or not the Lagrangian acceleration
  integer,save :: ilagr_vort ! compute or not the Lagrangian vorticity
  integer,save :: itraj_test ! save the first particle in ascci file traj_test
  real(kind=8),save :: tsave_part ! Time between outpout of entire fields.  
  character(len=100),save :: inicond_lagrangian ! initial positions of the particles
  character(len=100),save :: lagr_interp ! interpolation method

contains
   !-----------------------------------------------------------------
    include "init_lagrangian.f90"
    include "lagrangian_acc.f90"
    include "lagrangian_backup.f90"
    include "lagrangian_migration.f90"
    include "params_lagrangian.f90"
    include "lagrangian.f90"
   !---------------------------------------------------------------------------
    real(kind=8) function per_x(x,xl,dx)
      implicit none
      real(kind=8), intent (in) ::x, xl, dx
      real(kind=8) :: tmp
      tmp=x
      !if (xi(0)<0.d0) then
      !    if (tmp>xl/2.)  tmp = tmp-xl/2.
      !    if (tmp<-xl/2.) tmp = tmp+xl/2.
      !else
          if (tmp>xl)      tmp = tmp-xl
          if (tmp<0.d0)    tmp = tmp+xl        
      !endif
      per_x=tmp
      return
    end function per_x    
    !---------------------------------------------------------------------------

end module particles
