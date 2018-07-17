module module_insects_flusi
  use vars, only : xl, yl, zl, abort, &
  periodize_coordinate, cross, deg2rad, pi, rad2deg, root, norm2, &
  x0, y0, z0, startup_conditioner, rand_nbr, &
  itdrag, Integrals, nu

  ! interp2_nonper: we need this to interpolate wing thickness and corrugation
  use module_helpers, only : fseries_eval, hermite_eval, mpisum, interp2_nonper
  ! we need this to read from ini files (e.g. the wing kinematics or shape are read this way)
  use module_ini_files_parser_mpi

end module module_insects_flusi
