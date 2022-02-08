
subroutine post_force(help)
    use vars
    use p3dfft_wrapper
    use basic_operators
    use module_helpers
    use module_ini_files_parser_mpi
    use penalization ! mask, and us array
    use hdf5_wrapper

    implicit none
    logical, intent(in) :: help
    character(len=strlen) :: fname_mask, outfile, dummy
    character(len=strlen), allocatable :: fname_us(:), fname_u(:)
    real(kind=pr), dimension(:,:,:,:), allocatable :: u
    real(kind=pr) :: time, tmp(1)
    real(kind=pr), allocatable, DIMENSION(:) :: force
    integer :: ncomponents, j, i

    if (help.and.root) then
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "./flusi -p --force --components 3 --mask mask_0000.h5 --us usx_000.h5 &
        &[usy_00.h5] --u ux_00.h5 [uy_00.h5] -o forces.txt"
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "Parallel: Yes"
        return
    endif


    ! fetch parameters from command line call
    do i = 1, COMMAND_ARGUMENT_COUNT()
        call get_command_argument(i,dummy)

        select case (dummy)
        case ("--mask")
            call get_command_argument(i+1, fname_mask)
            if (root) write(*,*) "mask function read from=", fname_mask
            call check_file_exists(fname_mask)

        case ("--us")
            do j = 1, ncomponents
                call get_command_argument(i+j, fname_us(j))
                call check_file_exists(fname_us(j))
            enddo
            if (root) write(*,*) "solid velocity fields=", fname_us

        case ("--u")
            do j = 1, ncomponents
                call get_command_argument(i+j, fname_u(j))
                call check_file_exists(fname_u(j))
            enddo
            if (root) write(*,*) "velocity fields=", fname_u

        case ("-o")
            call get_command_argument(i+1,outfile)
            if (root) write(*,*) "will write output to=", outfile

        case ("--components")
            call get_command_argument(i+1,dummy)
            read(dummy,*) ncomponents

            if (root) write(*,*) "Expect ncomponents=", ncomponents

            allocate(fname_us(1:ncomponents))
            allocate(fname_u(1:ncomponents))
            allocate(force(1:ncomponents))

        end select
    enddo

    ! get file size etc
    call fetch_attributes( fname_mask, nx, ny, nz, xl, yl, zl, time, nu, origin )
    ! fetch value of penalization parameter from mask file. It should be stored
    ! there for every version of flusi. Note we do not need to divide by eps here, since
    ! this is done in RHS. It is the first time we need
    ! the epsi parameter since 2014 (now is 2018), so it was not included in
    ! fetch_attributes
    call read_attribute( fname_mask, "mask", "epsi", tmp)
    eps = tmp(1)
    if (root) write(*,*) "penalization parameter is eps=", eps
    ! initialization parallel module (no FFTS)
    call decomposition_initialize()

    allocate( u( ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3 ) )
    allocate( us( ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3 ) )
    allocate( mask( ra(1):rb(1),ra(2):rb(2),ra(3):rb(3) ) )
    ! allocate( mask_color( ra(1):rb(1),ra(2):rb(2),ra(3):rb(3) ) )

    call read_single_file( fname_mask, mask )
    do j = 1, ncomponents
        call read_single_file( fname_u(j), u(:,:,:,j) )
        call read_single_file( fname_us(j), us(:,:,:,j) )
    enddo

    mask = mask / eps

    do j = 1, ncomponents
        force(j) = sum( mask*(u(:,:,:,j)-us(:,:,:,j)) )
    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE, force, ncomponents, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, i)

    force = force * dx*dy*dz

    if (root) then
        write(*,*) 'writing temporal coefficients ', outfile
        open(14, file=outfile, status='replace')
        write(14,'(10(es15.8,1x))') time, force
        close(14)
    end if
end subroutine
