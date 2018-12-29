!-------------------------------------------------------------------------------
! ./flusi --postprocess --time-avg file_list.txt avgx_0000.h5
! Reads in a list of files from a file, then loads one file after the other and
! computes the average field, which is then stored in the specified file.
!-------------------------------------------------------------------------------
subroutine POD(help)
    use vars
    use p3dfft_wrapper
    use basic_operators
    use module_helpers
    use module_ini_files_parser_mpi
    implicit none
    logical, intent(in) :: help
    character(len=strlen) :: fnamex_list, fname_this, fnamey_list, fnamez_list, dummy
    character(len=strlen) :: fnamex, fnamey, fnamez
    real(kind=pr), dimension(:,:,:), allocatable :: field_avg, field
    real(kind=pr), dimension(:,:), allocatable :: X_data, POD_modes, a_coefs
    real(kind=pr), dimension(:), allocatable :: X_mean
    DOUBLE PRECISION, dimension(:,:), allocatable :: C, V, D
    DOUBLE PRECISION, dimension(:), allocatable :: eigenvalues, work
    integer :: ix, iy ,iz, io_error, i,j, N_modes, N_snapshots, info, it, dim, npoints
    real(kind=pr) :: time, a, norm
    LOGICAL :: vector

    if (help.and.root) then
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "./flusi -p --POD [--vector] [--2D | --3D] --list file_list.txt [list_uy.txt] [list_uz.txt] --modes 10"
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) " --list: a TXT file which simply contains the list of snapshots to read, one file per line"
        write(*,*) " "
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "Parallel: Yes"
        return
    endif

    ! defaults:
    vector = .false.
    dim = 2

    ! fetch parameters from command line call
    do i = 1, COMMAND_ARGUMENT_COUNT()
        call get_command_argument(i,dummy)
        select case (dummy)
        case ("--list")
            if (vector) then
                ! vector
                call get_command_argument(i+1,fnamex_list)
                call get_command_argument(i+2,fnamey_list)
                call get_command_argument(i+3,fnamez_list)
            else
                ! scalar
                call get_command_argument(i+1,fnamex_list)
            endif

        case ("--modes")
            call get_command_argument(i+1,dummy)
            read(dummy,*) N_modes
            if (root) write(*,*) "Will save N_modes=", N_modes

        case ("--vector")
            vector = .true.

        case ("--scalar")
            vector = .false.

        case ("--2D")
            dim = 2

        case ("--3D")
            dim = 3

        end select
    enddo

    !-----------------------------------------------------------------------------
    ! check if input file exists, the file contains the list of h5 files to be avg
    !-----------------------------------------------------------------------------
    call check_file_exists ( fnamex_list )
    if (vector) call check_file_exists ( fnamey_list )
    if (vector .and. dim==3) call check_file_exists ( fnamez_list )

    if (root) write(*,*) "Reading list of files from "//fnamex_list
    if (root) write(*,*) fnamex_list, fnamey_list, fnamez_list, vector, dim, N_modes

    !-----------------------------------------------------------------------------
    ! read in the file, loop over lines
    !-----------------------------------------------------------------------------
    call count_lines_in_ascii_file_mpi(fnamex_list, N_snapshots, n_header=0)

    if (root) write(*,*) "Reading N_snapshots=", N_snapshots


    allocate( C(1:N_snapshots,1:N_snapshots) )
    allocate( D(1:N_snapshots,1:N_snapshots) )
    allocate( V(1:N_snapshots,1:N_snapshots) )
    allocate( eigenvalues(1:N_snapshots) )
    allocate( a_coefs(1:N_snapshots,1:N_snapshots) )
    allocate( work(1:5*N_snapshots) )


    open( unit=14, file=fnamex_list, action='read', status='old' )
    if (vector) then
        open( unit=15, file=fnamey_list, action='read', status='old' )
    endif
    if (vector .and. dim==3) then
        open( unit=16, file=fnamez_list, action='read', status='old' )
    endif

    io_error = 0
    i = 1
    do while (i <= N_snapshots)
        read (14,'(A)', iostat=io_error) fnamex
        call check_file_exists ( fnamex )

        if (vector) then
            read (15,'(A)', iostat=io_error) fnamey
            call check_file_exists ( fnamey )
        endif

        if (vector .and. dim==3) then
            read (16,'(A)', iostat=io_error) fnamez
            call check_file_exists ( fnamez )
        endif

        ! initialization is done after first read.
        if ( i == 1 ) then
            ! get file size etc
            call fetch_attributes( fnamex, nx, ny, nz, xl, yl, zl, time, nu, origin )
            ! initialization parallel module (no FFTS)
            call decomposition_initialize()

            ! size of a flattened snapshot
            npoints = (rb(1)-ra(1)+1) * (rb(2)-ra(2)+1) * (rb(3)-ra(3)+1)

            if (decomposition /= "1D") call abort(28122018,"I think this module works only for 1D MPI decomp (scalar products)")

            ! memory
            allocate( field( ra(1):rb(1),ra(2):rb(2),ra(3):rb(3) ) )

            if (vector) then
                allocate( X_data( npoints * dim, 1:N_snapshots) )
                allocate( POD_modes( npoints * dim, 1:N_snapshots) )
                allocate( X_mean( npoints * dim ) )

            else
                allocate( X_data( npoints, 1:N_snapshots) )
                allocate( POD_modes( npoints, 1:N_snapshots) )
                allocate( X_mean( npoints ) )

            endif

        endif

        ! read the field from file
        call read_single_file( fnamex, field )

        ! add it to the data array
        X_data(1:npoints,i) = reshape(field, (/npoints/) )

        if (vector) then
            ! read the field from file
            call read_single_file( fnamey, field )

            ! add it to the data array
            X_data(npoints+1:2*npoints,i) = reshape(field, (/npoints/) )
        endif

        if (vector .and. dim==3) then
            ! read the field from file
            call read_single_file( fnamez, field )

            ! add it to the data array
            X_data(2*npoints+1:3*npoints,i) = reshape(field, (/npoints/) )
        endif

        if (root) write(*,*) "filled snapshot slot", i
        i = i+1
    enddo

    close (14)
    if (vector) close (15)
    if (vector .and. dim==3) close (16)

    !---------------------------------------------------------------------------
    ! compute fluctuations (remove mean)
    !---------------------------------------------------------------------------

    ! compute mean of data
    X_mean = 0.0_pr
    do i = 1, N_snapshots
        X_mean = X_mean + X_data(:,i)
    enddo
    X_mean = X_mean / dble(N_snapshots)

    ! POD acts on fluctuations, so remove the mean.
    do i = 1, N_snapshots
        X_data(:,i) = X_data(:,i) - X_mean
    enddo

    ! divide by number of snapshots (eqn. 3.29 on p. 11 of Luchtenberg, Noak.)
    X_data = X_data / dble(N_snapshots)

    !---------------------------------------------------------------------------
    ! construction of covariance matrix from snapsot data
    !---------------------------------------------------------------------------
    ! compute matrix C
    C = matmul( transpose(X_data), X_data )
    ! do i = 1, N_snapshots
    !     do j = 1, N_snapshots
    !         C(i,j) = sum( X_data(:,i)*X_data(:,j) )
    !     enddo
    ! enddo

    ! normalization of C Matrix
    C = C / dble( nx*ny*nz )
    call MPI_ALLREDUCE(MPI_IN_PLACE, C, N_snapshots**2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)


    ! if (root) then
    !     open(14, file='C_matrix_fortran.txt', status='replace')
    !     do i = 1, N_snapshots
    !         write(14,'(400(es15.8,1x))') C(i,:)
    !     enddo
    !     close(14)
    ! end if

    !---------------------------------------------------------------------------
    ! eigenvalues of covariance matrix
    !---------------------------------------------------------------------------
    call DSYEV('V', 'U', N_snapshots, C, N_snapshots, eigenvalues, work, 5*N_snapshots, info)
    ! as in matlab the eigenvalues are sorted in ascending order...
    ! on output, C now contains the eigenvectors:
    V = C

    if (root) write(*,*) "info=", info
    if (info /= 0) call abort(333,"The eigenvalue solver failed...")


    if (root) then
        write(*,*) "----v eigenvalues v-----"
        do i = 1, N_snapshots
            write(*,'(1(es12.4,1x))') eigenvalues(i)
        enddo
        write(*,*) "----^ eigenvalues ^-----"
    endif

    !---------------------------------------------------------------------------
    ! construct POD basis functions (modes)
    !---------------------------------------------------------------------------
    ! eqn. 31 from AIAA review

    POD_modes = matmul(X_data, V)

    ! do i = 1, size(POD_modes, 1)
    !     do j = 1, N_snapshots
    !         POD_modes(i,j) = sum( X_data(i,:)*V(:,j) )
    !     enddo
    ! enddo

    do i = 1, N_snapshots
        ! ! note the division by sqrt(lambda) might cause numerical instabilities
        ! ! if the eigenvalue is close to zero
        ! if ( eigenvalues(i) > 1.0e-9_pr ) then
        !     POD_modes(:,i) = POD_modes(:,i) / sqrt(eigenvalues(i))
        ! endif

        ! for some reason the modes are not normalized here, so we take care of that
        ! I think the reason is the scalar product, which needs to be scaled by nx*ny*nz (and NOT npoints)
        norm = sqrt( sum(POD_modes(:,i)**2) )
        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)

        if (norm > 1.0e-8_pr) then
            POD_modes(:,i) = POD_modes(:,i) / norm
        endif

    enddo

    !---------------------------------------------------------------------------
    ! save modes to disk
    !---------------------------------------------------------------------------
    do i = 1, N_modes
        field = reshape( POD_modes(1:npoints,N_snapshots-i+1), &
        (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

        ! create filename
        write(fname_this,'("modex_",i3.3,".h5")') i

        call save_field_hdf5 ( dble(i), fname_this, field )
    enddo

    if (vector) then
        do i = 1, N_modes
            field = reshape( POD_modes(npoints+1:2*npoints,N_snapshots-i+1), &
            (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

            ! create filename
            write(fname_this,'("modey_",i3.3,".h5")') i

            call save_field_hdf5 ( dble(i), fname_this, field )
        enddo
    endif

    if (vector .and. dim==3) then
        do i = 1, N_modes
            field = reshape( POD_modes(2*npoints+1:3*npoints,N_snapshots-i+1), &
            (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

            ! create filename
            write(fname_this,'("modez_",i3.3,".h5")') i

            call save_field_hdf5 ( dble(i), fname_this, field )
        enddo
    endif


    !---------------------------------------------------------------------------
    ! temporal coefficients
    !---------------------------------------------------------------------------
    a_coefs = 0.0_pr
    do it = 1, N_snapshots
        do i = 1, N_modes
            ! scalar product (inner product)
            a_coefs(it,i) = sum( X_data(:, it) * POD_modes(:, N_snapshots-i+1) )
        enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE, a_coefs, N_snapshots**2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)

    if (root) then
        write(*,*) 'writing temporal coefficients a_coefs.txt'
        open(14, file='a_coefs.txt', status='replace')
        do i = 1, N_snapshots
            write(14,'(400(es15.8,1x))') a_coefs(i,1:N_modes)
        enddo
        close(14)
    end if

    !---------------------------------------------------------------------------
    ! reconstruction using N_modes
    !---------------------------------------------------------------------------
    it = 1

    ! ------------------ x-component ---------------------
    field = 0.0_pr
    do i = 1, N_modes
        field = field + a_coefs(it,i) * reshape( POD_modes(1:npoints, N_snapshots-i+1), &
        (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )
    enddo

    ! create filename
    write(fname_this,'("reconstructionx_",i3.3,".h5")') N_modes
    call save_field_hdf5 ( dble(it), fname_this, field )

    ! for comparison, also save original data (which is of course only the fluctuating
    ! part of the solution)
    field = reshape( X_data(1:npoints, it), (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

    write(fname_this,'("originalx_",i3.3,".h5")') N_modes
    call save_field_hdf5 ( dble(it), fname_this, field )


    ! ------------------ y-component ---------------------
    if (vector) then
        field = 0.0_pr
        do i = 1, N_modes
            field = field + a_coefs(it,i) * reshape( POD_modes(npoints+1:2*npoints, N_snapshots-i+1), &
            (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )
        enddo

        ! create filename
        write(fname_this,'("reconstructiony_",i3.3,".h5")') N_modes
        call save_field_hdf5 ( dble(it), fname_this, field )

        ! for comparison, also save original data (which is of course only the fluctuating
        ! part of the solution)
        field = reshape( X_data(npoints+1:2*npoints, it), (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

        write(fname_this,'("originaly_",i3.3,".h5")') N_modes
        call save_field_hdf5 ( dble(it), fname_this, field )
    endif


    ! ------------------ z-component ---------------------
    if (vector .and. dim==3) then
        field = 0.0_pr
        do i = 1, N_modes
            field = field + a_coefs(it,i) * reshape( POD_modes(2*npoints+1:3*npoints, N_snapshots-i+1), &
            (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )
        enddo

        ! create filename
        write(fname_this,'("reconstructionz_",i3.3,".h5")') N_modes
        call save_field_hdf5 ( dble(it), fname_this, field )

        ! for comparison, also save original data (which is of course only the fluctuating
        ! part of the solution)
        field = reshape( X_data(2*npoints+1:3*npoints, it), (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

        write(fname_this,'("originalz_",i3.3,".h5")') N_modes
        call save_field_hdf5 ( dble(it), fname_this, field )
    endif


end subroutine POD
