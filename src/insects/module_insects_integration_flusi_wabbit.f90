module module_insects_integration_flusi_wabbit
  use vars, only : abort, periodize_coordinate, cross, deg2rad, pi, rad2deg, norm2, &
  startup_conditioner, rand_nbr, Integrals, x0, y0, z0, rk => pr

  ! interp2_nonper: we need this to interpolate wing thickness and corrugation
  use module_helpers
  ! we need this to read from ini files (e.g. the wing kinematics or shape are read this way)
  use module_ini_files_parser_mpi

implicit none

  logical, parameter :: grid_time_dependent = .false.

end module module_insects_integration_flusi_wabbit
