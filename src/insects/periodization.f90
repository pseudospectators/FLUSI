! given a point x, check if it lies in the computational domain centered at zero
! (note: we assume [-xl/2...+xl/2] size this is useful for insects )
function periodize_coordinate(x_glob)
  real(kind=pr) :: x_glob(1:3)
  real(kind=pr),dimension(1:3) :: periodize_coordinate
  periodize_coordinate = x_glob

  if (periodic) then
    if (x_glob(1)<-xl/2.0) periodize_coordinate(1)=x_glob(1)+xl
    if (x_glob(2)<-yl/2.0) periodize_coordinate(2)=x_glob(2)+yl
    if (x_glob(3)<-zl/2.0) periodize_coordinate(3)=x_glob(3)+zl

    if (x_glob(1)>xl/2.0) periodize_coordinate(1)=x_glob(1)-xl
    if (x_glob(2)>yl/2.0) periodize_coordinate(2)=x_glob(2)-yl
    if (x_glob(3)>zl/2.0) periodize_coordinate(3)=x_glob(3)-zl
  else
    periodize_coordinate = x_glob
  endif

end function
