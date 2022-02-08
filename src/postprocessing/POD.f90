subroutine POD(help)
    use vars
    use p3dfft_wrapper
    use basic_operators
    use module_helpers
    use module_ini_files_parser_mpi
    implicit none
    logical, intent(in) :: help
    character(len=strlen) :: fnamex_list, fname_this, fnamey_list, fnamez_list, dummy
    character(len=strlen) :: fnamex, fnamey, fnamez, timesteps_list, fname_reconst_base
    character(len=strlen) :: fname_list(1:3), fname_acoefs
    real(kind=pr), dimension(:,:,:), allocatable :: field_avg, field
    real(kind=pr), dimension(:,:), allocatable :: X_data, POD_modes, a_coefs
    real(kind=pr), dimension(:), allocatable :: X_mean, times
    DOUBLE PRECISION, dimension(:,:), allocatable :: C, V, D
    DOUBLE PRECISION, dimension(:), allocatable :: eigenvalues, work
    integer :: ix, iy ,iz, io_error, i,j, N_modes, N_snapshots, info, it, dim, npoints, k
    integer :: ncomponents, nmodes_to_read, nmodes_file
    real(kind=pr) :: time, a, norm
    LOGICAL :: vector, reconstruct_all_time_steps, reconstruct_list, only_recon, total_not_fluctuations, only_modes
    logical :: save_original

    if (help.and.root) then
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "./flusi -p --POD --components 3 --list file_list.txt [list_uy.txt] [list_uz.txt] "
        write(*,*) "--modes 10 --reconstruct-all-time-steps --total-not-fluctuations"
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) " Snapshot POD of flow data (vorticity, velocity, etc.) "
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) " --list: a TXT file which simply contains the list of snapshots to read, one file per line"
        write(*,*) " --reconstruct-all-time-steps"
        write(*,*) " --reconstruct-list timesteps.txt"
        write(*,*) " --only-modes"
        write(*,*) " --only-reconstruction a_coefs.txt (reconstruction from 34 modes, read from HDD)"
        write(*,*) " --modes 10"
        write(*,*) " --save-original"
        write(*,*) "    if specified, code also stores original fluctuations (e.g. input data) at reconstruction"
        write(*,*) "    time steps. In the case of --total-not-fluctuations, this would be useless"
        write(*,*) " --total-not-fluctuations"
        write(*,*) "    do not remove mean from data if this flag is set. in the literature, POD always removes the"
        write(*,*) "    ensemble average."
        write(*,*) " --reconstruction-name vel-fluct"
        write(*,*) "    If the code is reconstructing the field from POD modes (i.e. not running --only-modes)"
        write(*,*) "    then you can choose the filename for reconstructions here. Naming scheme is"
        write(*,*) "    [NAME]-[MODES]-[COMPONENT]_[TIMESTEP].h5"
        write(*,*) "    so e.g. vel-fluct-30-1_000.h5 if name='vel-fluct' "
        write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        write(*,*) "Parallel: Yes"
        return
    endif

    ! defaults:
    ncomponents = 1
    reconstruct_all_time_steps = .false.
    reconstruct_list = .false.
    only_recon = .false.
    total_not_fluctuations = .false.
    only_modes = .false.
    save_original = .false.
    fname_reconst_base = "reconstruction"

    ! fetch parameters from command line call
    do i = 1, COMMAND_ARGUMENT_COUNT()
        call get_command_argument(i,dummy)
        select case (dummy)
        case ("--save-original")
            save_original = .true.
            if (root) write(*,*) "We do save original data."

        case ("--total-not-fluctuations")
            total_not_fluctuations = .true.
            if (root) write(*,*) "We do NOT remove the mean and work on TOTAL field"

        case ("--reconstruction-name")
            call get_command_argument(i+1,fname_reconst_base)
            if (root) write(*,*) "Basename for reconstruction is", fname_reconst_base

        case ("--only-modes")
            only_modes = .true.
            if (root) write(*,*) "We compute only modes and do NOT reconstruct"

        case ("--only-reconstruction")
            only_recon = .true.

            call get_command_argument(i+1, fname_acoefs)
            if (root) write(*,*) "reconstruction from acoefs=", fname_acoefs

        case ("--list")
            do j = 1, ncomponents
                call get_command_argument(i+j, fname_list(j))
            enddo

        case ("--reconstruct-all-time-steps")
            reconstruct_all_time_steps = .true.

        case ("--modes")
            call get_command_argument(i+1,dummy)
            read(dummy,*) N_modes
            if (root) write(*,*) "Will save N_modes=", N_modes

        case ("--reconstruct-list")
            call get_command_argument(i+1,timesteps_list)
            reconstruct_list = .true.
            reconstruct_all_time_steps = .false.
            if (root) write(*,*) "Will reconstruct at time steps given by ", timesteps_list
            call check_file_exists(timesteps_list)

        case ("--components")
            call get_command_argument(i+1,dummy)
            read(dummy,*) ncomponents
            if (root) write(*,*) "Expect ncomponents=", ncomponents

        end select
    enddo

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! VARIANT A: read data, compute POD
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (.not. only_recon) then
        !-----------------------------------------------------------------------------
        ! check if input file exists, the file contains the list of h5 files to be avg
        !-----------------------------------------------------------------------------
        do j = 1, ncomponents
            call check_file_exists ( fname_list(j) )
            if (root) write(*,*) "Reading list of files from "//fname_list(j)
        enddo


        !-----------------------------------------------------------------------------
        ! read in the file, loop over lines
        !-----------------------------------------------------------------------------
        call count_lines_in_ascii_file_mpi(fname_list(1), N_snapshots, n_header=0)

        if (root) write(*,*) "Reading N_snapshots=", N_snapshots

        allocate( C(1:N_snapshots,1:N_snapshots) )
        allocate( D(1:N_snapshots,1:N_snapshots) )
        allocate( V(1:N_snapshots,1:N_snapshots) )
        allocate( eigenvalues(1:N_snapshots) )
        allocate( a_coefs(1:N_snapshots,1:N_snapshots) )
        allocate( work(1:5*N_snapshots) )
        allocate( times(1:N_snapshots) )

        do j = 1, ncomponents
            open( unit=10+j, file=fname_list(j), action='read', status='old' )
        enddo

        io_error = 0
        i = 1
        do while (i <= N_snapshots)
            do j = 1, ncomponents
                read (10+j, '(A)', iostat=io_error) fnamex
                call check_file_exists ( fnamex )

                ! initialization is done after first read.
                if ( i == 1 .and. j == 1 ) then
                    ! get file size etc
                    call fetch_attributes( fnamex, nx, ny, nz, xl, yl, zl, time, nu, origin )
                    ! initialization parallel module (no FFTS)
                    call decomposition_initialize()

                    ! size of a flattened snapshot
                    npoints = (rb(1)-ra(1)+1) * (rb(2)-ra(2)+1) * (rb(3)-ra(3)+1)

                    if (decomposition /= "1D") then
                        call abort(28122018,"I think this module works only for 1D MPI decomp (scalar products)")
                    endif

                    ! memory for one field
                    allocate( field( ra(1):rb(1),ra(2):rb(2),ra(3):rb(3) ) )

                    allocate( X_data( 1:npoints * ncomponents, 1:N_snapshots) )
                    allocate( POD_modes( 1:npoints * ncomponents, 1:N_snapshots) )
                    allocate( X_mean( 1:npoints * ncomponents ) )
                endif

                ! read the field from file
                call read_single_file( fnamex, field )

                ! add it to the data array
                X_data( 1+(j-1)*npoints:(j)*npoints, i) = reshape(field, (/npoints/) )
            enddo
            if (root) write(*,*) "filled snapshot slot", i
            i = i+1
        enddo

        do j = 1, ncomponents
            close(10+j)
        enddo

        !---------------------------------------------------------------------------
        ! compute fluctuations (remove mean)
        !---------------------------------------------------------------------------
        if ( .not. total_not_fluctuations ) then
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
        endif

        !---------------------------------------------------------------------------
        ! construction of covariance matrix from snapshot data
        !---------------------------------------------------------------------------
        ! compute matrix C
        C = matmul( transpose(X_data), X_data )

        ! divide by number of snapshots (eqn. 3.29 on p. 11 of Luchtenberg, Noak.)
        C = C / dble(N_snapshots)

        call MPI_ALLREDUCE(MPI_IN_PLACE, C, N_snapshots**2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)

        ! normalization of C Matrix
        C = C / dble( nx*ny*nz )


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

        if (root) write(*,*) "Constructing POD modes (X*V)"
        POD_modes = matmul(X_data, V)


        if (root) write(*,*) "Constructing POD modes (sqrt(lambda))"
        do i = 1, N_snapshots
            ! note the division by sqrt(lambda) might cause numerical instabilities
            ! if the eigenvalue is close to zero
            if ( eigenvalues(i) > 1.0e-9_pr ) then
                POD_modes(:,i) = POD_modes(:,i) / sqrt( dble(N_snapshots)*eigenvalues(i))
                ! POD_modes(:,i) = POD_modes(:,i) / sqrt(eigenvalues(i))
            endif

            ! ! for some reason the modes are not normalized here, so we take care of that
            ! ! I think the reason is the scalar product, which needs to be scaled by nx*ny*nz (and NOT npoints)
            ! norm = sqrt( sum(POD_modes(:,i)**2) / dble(nx*ny*nz) )
            !
            ! call MPI_ALLREDUCE(MPI_IN_PLACE, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
            !
            ! if (root) write(*,*) i, norm
            !
            ! if (norm > 1.0e-8_pr) then
            !     POD_modes(:,i) = POD_modes(:,i) / norm
            ! endif
        enddo

        !---------------------------------------------------------------------------
        ! save modes to disk
        !---------------------------------------------------------------------------
        if (root) write(*,*) "Saving POD modes to disk"

        do i = 1, N_modes
            do j = 1, ncomponents
                field = reshape( POD_modes(1+(j-1)*npoints:(j)*npoints, N_snapshots-i+1), &
                (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

                ! create filename
                write(fname_this, '("mode",i1,"_",i3.3,".h5")') j, i

                call save_field_hdf5 ( dble(i), fname_this, field )
            enddo
        enddo

        !---------------------------------------------------------------------------
        ! temporal coefficients
        !---------------------------------------------------------------------------
        if (root) write(*,*) "Computing temporal coefficients a"

        a_coefs = 0.0_pr
        do it = 1, N_snapshots
            do i = 1, N_modes
                ! scalar product (inner product)
                a_coefs(it,i) = sum( X_data(:, it) * POD_modes(:, N_snapshots-i+1) )
            enddo
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, a_coefs, N_snapshots**2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
        a_coefs = a_coefs / dble(nx*ny*nz)

        if (root) then
            write(*,*) 'writing temporal coefficients a_coefs.txt'
            open(14, file='a_coefs.txt', status='replace')
            do i = 1, N_snapshots
                write(14,'(400(es15.8,1x))') a_coefs(i, 1:N_modes)
            enddo
            close(14)
        end if


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! VARIANT B: read previously computed POD from files
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else

        !**************** read-from-file-mode ****************

        call count_lines_in_ascii_file_mpi(fname_acoefs, N_snapshots, n_header=0)
        call count_cols_in_ascii_file_mpi(fname_acoefs, nmodes_file, n_header=0)

        if (root) write(*,*) "a_coefs.txt:", N_snapshots, nmodes_file

        allocate( a_coefs(1:N_snapshots, 1:nmodes_file))
        call read_array_from_ascii_file_mpi(fname_acoefs, a_coefs, n_header=0)

        if (root) then
            write(*,*) "~~~~v      ", fname_acoefs
            do k = 1, N_snapshots
                write(*,'(400(es15.8,1x))') a_coefs(k, :)
            enddo
            write(*,*) "~~~~^      ", fname_acoefs
        endif

        if (N_snapshots < N_modes) then
            call abort(7771,"You want to reconstruct using more modes than the POD originally used.")
        endif

        ! read modes, but only as many as we require for reconstruction
        do i = 1, N_modes
            do j = 1, ncomponents
                ! create filename
                write(fname_this, '("mode",i1,"_",i3.3,".h5")') j, i

                if (root) write(*,*) "reading mode: ", fname_this

                if (i==1 .and. j==1) then
                    ! get file size etc
                    call fetch_attributes( fname_this, nx, ny, nz, xl, yl, zl, time, nu, origin )
                    ! initialization parallel module (no FFTS)
                    call decomposition_initialize()

                    ! size of a flattened snapshot
                    npoints = (rb(1)-ra(1)+1) * (rb(2)-ra(2)+1) * (rb(3)-ra(3)+1)

                    if (decomposition /= "1D") then
                        call abort(28122018,"I think this module works only for 1D MPI decomp (scalar products)")
                    endif

                    ! memory for one field
                    allocate( field( ra(1):rb(1),ra(2):rb(2),ra(3):rb(3) ) )
                    allocate( POD_modes( 1:npoints * ncomponents, 1:N_snapshots) )
                endif

                ! read the mode
                call read_single_file( fname_this, field )

                ! sort it in the array (note this step is not necessary but simplifies coding)
                POD_modes( 1+(j-1)*npoints:(j)*npoints, N_snapshots-i+1) = reshape(field, (/npoints/) )

            enddo
        enddo

    endif


    ! at this point, we have the POD_modes and their temporal coefficients ready

    ! if not reconstructing, we're done now
    if (only_modes) return


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! RECONSTRUCTION
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !---------------------------------------------------------------------------
    ! reconstruction using N_modes
    !---------------------------------------------------------------------------
    ! which time step to reconstruct at
    if (reconstruct_all_time_steps .and. .not. reconstruct_list) then
        do it = 1, N_snapshots
            do j = 1, ncomponents
                field = 0.0_pr

                do i = 1, N_modes
                    field = field + a_coefs(it,i) * reshape( POD_modes(1+(j-1)*npoints:(j)*npoints, N_snapshots-i+1), &
                    (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )
                enddo

                ! create filename
                write(fname_this,'(A,"-",i2.2,"-",i1,"_",i3.3,".h5")') trim(adjustl(fname_reconst_base)), N_modes, j, it
                call save_field_hdf5 ( dble(it), fname_this, field )

                if (.not. only_recon .and. save_original) then
                    ! for comparison, also save original data (which is of course only the fluctuating
                    ! part of the solution). can only be done when computing POD, not on reconstruction
                    field = reshape( X_data(1+(j-1)*npoints:(j)*npoints, it), &
                    (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

                    write(fname_this,'("original",i2.2,"-",i1,"_",i3.3,".h5")') N_modes, j, it
                    call save_field_hdf5 ( dble(it), fname_this, field )
                endif

            enddo
        enddo

    !---------------------------------------------------------------------------
    ! reconstruction of a list of time steps
    !---------------------------------------------------------------------------
    elseif ( .not. reconstruct_all_time_steps .and. reconstruct_list ) then

        io_error = 0
        open(17, file=timesteps_list, action='read', status='old' )

        do while (io_error == 0)
            read(17, *, iostat=io_error) it

            if (io_error == 0) then

                if (root) write(*,*) "reconstruction it=", it, "using", N_modes, ncomponents

                if (it<1 .or. it>size(a_coefs,1)) then
                    write(*,*) "it=", it, size(a_coefs,1)
                    call abort(8881,"You request reconstruction at an invalid time step.")
                endif

                do j = 1, ncomponents
                    field = 0.0_pr

                    do i = 1, N_modes
                        field = field + a_coefs(it,i) * reshape( POD_modes(1+(j-1)*npoints:(j)*npoints, N_snapshots-i+1), &
                        (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )
                    enddo

                    ! create filename
                    write(fname_this,'(A,"-",i2.2,"-",i1,"_",i3.3,".h5")') trim(adjustl(fname_reconst_base)), N_modes, j, it
                    call save_field_hdf5 ( dble(it), fname_this, field )

                    if (.not. only_recon .and. save_original) then
                        ! for comparison, also save original data (which is of course only the fluctuating
                        ! part of the solution). can only be done when computing POD, not on reconstruction
                        field = reshape( X_data(1+(j-1)*npoints:(j)*npoints, it), &
                        (/(rb(1)-ra(1)+1), (rb(2)-ra(2)+1), (rb(3)-ra(3)+1)/) )

                        write(fname_this,'("original",i2.2,"-",i1,"_",i3.3,".h5")') N_modes, j, it
                        call save_field_hdf5 ( dble(it), fname_this, field )
                    endif
                enddo
            endif
        enddo

        close(17)
    endif
end subroutine POD
