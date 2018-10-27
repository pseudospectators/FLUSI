module krylov_module
    use vars

contains
    subroutine krylov(time,it,dt0,dt,u,uk,nlk,vort,work,workc,expvis,press,&
        scalars,scalars_rhs,Insect,beams)
        use mpi
        use vars
        use p3dfft_wrapper
        use solid_model
        use module_insects
        use basic_operators
        implicit none

        real(kind=pr),intent(inout) :: time,dt,dt0
        integer,intent(in) :: it
        complex(kind=pr),intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
        complex(kind=pr),intent(inout) :: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:nrhs-1)
        complex(kind=pr),intent(inout) :: workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
        real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
        real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
        real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
        real(kind=pr),intent(inout) :: expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
        ! pressure array. this is with ghost points for interpolation
        real(kind=pr),intent(inout) :: press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
        real(kind=pr),intent(inout) :: scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
        real(kind=pr),intent(inout) :: scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
        type(solid),dimension(1:nbeams),intent(inout) :: beams
        type(diptera),intent(inout) :: Insect

        real(kind=pr), allocatable, save :: H(:,:), phiMat(:,:), H_tmp(:,:), rhs(:,:,:,:,:), uold(:,:,:,:)
        complex(kind=pr), allocatable, save :: uktmp(:,:,:,:)
        integer :: M_iter
        integer :: i, j, k, l, iter
        real(kind=pr) :: normv, eps2, beta
        real(kind=pr) :: h_klein, err, t0


        ! allocate matrices with largest admissible
        if (.not. allocated(H)) then
            allocate( H(M_max+2,M_max+2) )
            allocate( H_tmp(M_max+2,M_max+2) )
            allocate( phiMat(M_max+2,M_max+2) )
        endif

        if (.not.allocated(rhs)) allocate( rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq,1:M_max+3) )
        if (.not.allocated(uold)) allocate( uold(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq) )
        if (.not.allocated(uktmp)) allocate( uktmp(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq) )

        phiMat = 0.0_pr
        H = 0.0_pr


        call ifft3( ink=uk, outx=uold )
        call adjust_dt(time, uold, dt)


        ! compute norm "normv" of input state vector
        call kry_norm( uold, normv )
        if (normv < epsilon(1.0_pr)) normv = 1.0_pr
        eps2 = normv * sqrt(epsilon(1.0_pr))


        ! the very last slot (M+3) is the "reference right hand side"
        call cal_nlk( time, it, nlk(:,:,:,:,0), uk, u, vort, work, workc, press, scalars, scalars_rhs, Insect, beams)
        call dealias( nlk(:,:,:,:,0) )
        call ifft3( ink=nlk(:,:,:,:,0), outx=rhs(:,:,:,:,M_max+3) )



        ! compute norm "beta", which is the norm of the reference RHS evaluation
        ! NSF: this guy is called normuu
        call kry_norm( rhs(:,:,:,:,M_max+3), beta )
        if (beta < epsilon(1.0_pr)) beta = 1.0_pr


        ! start iteration, fill first slot
        rhs(:,:,:,:,1) = rhs(:,:,:,:,M_max+3) / beta

        !**************************************!
        !*** begin interations              ***!
        !**************************************!
        ! we loop over the full size of subspace, then exit prematurely
        ! if possible.
        do M_iter = 1, M_max

            ! perturbed right hand side is first-to-last (M+2) slot
            rhs(:,:,:,:,M_max+2) = uold(:,:,:,:) + eps2 * rhs(:,:,:,:,M_iter)

            ! call RHS with perturbed state vector, stored in slot (M_max+1)
            call fft3( inx=rhs(:,:,:,:,M_max+2), outk=uktmp )
            call cal_nlk( time, it, nlk(:,:,:,:,0), uktmp, u, vort, work, workc, press, scalars, scalars_rhs, Insect, beams)
            call dealias( nlk(:,:,:,:,0) )
            call ifft3( ink=nlk(:,:,:,:,0), outx=rhs(:,:,:,:,M_max+1) )

            ! linearization of RHS slot (M_max+1)
            rhs(:,:,:,:,M_max+1) = ( rhs(:,:,:,:,M_max+1) - rhs(:,:,:,:,M_max+3) ) / eps2

            ! --- inner loop ----
            ! --- ARNOLDI ---
            do iter = 1, M_iter
                call scalarproduct( rhs(:,:,:,:,iter), rhs(:,:,:,:,M_max+1), H(iter, M_iter) )

                rhs(:,:,:,:,M_max+1) = rhs(:,:,:,:,M_max+1) - H(iter,M_iter) * rhs(:,:,:,:,iter)
            enddo
            ! end of inner i =1:j loop
            call kry_norm( rhs(:,:,:,:,M_max+1), H(M_iter+1,M_iter) )

            rhs(:,:,:,:,M_iter+1) = rhs(:,:,:,:,M_max+1) / H(M_iter+1,M_iter)

            ! if this is the last iteration, we compute the H matrix and the matrix exponential
            ! and use it to get an error estimate. If the error seems okay, we are done and can
            ! compute the new time step, otherwise we increase M by one and retry.
            ! if we use the dynamic method, we evaluate the error after every iteration to see if
            ! we're good to go.
            h_klein    = H(M_iter+1,M_iter)

            ! create a copy of the H matrix with the right dimension
            H_tmp      = 0.0_pr
            H_tmp(1:M_iter, 1:M_iter) = H(1:M_iter, 1:M_iter)
            H_tmp(M_iter+1,M_iter)    = 0.0_pr
            H_tmp(1,M_iter+1)         = 1.0_pr
            H_tmp(M_iter+1,M_iter+2)  = 1.0_pr


            ! compute matrix exponential
            phiMat = 0.0_pr
            phiMat(1:M_iter+2, 1:M_iter+2) = expM_pade( dt*H_tmp(1:M_iter+2, 1:M_iter+2) )
            phiMat(M_iter+1, M_iter+1)     = h_klein*phiMat(M_iter, M_iter+2)


            ! *** Error estimate ***!
            err = abs( beta*phiMat(M_iter+1,M_iter+1) )

            if ( M_iter == M_max .and. err > krylov_err_threshold) then
                if(mpirank==0) write(*,'("Krylov: did M_iter=",i2," got err=",es12.4," which is larger than your &
                &threshold. As no more iter possible, results may be inaccurate.")') M_iter, err
            endif

            if (err < krylov_err_threshold .or. M_iter == M_max) then
                if(mpirank==0) write(*,*) "exit krylov err=", err, "t=", time, dt, M_iter
                exit
            endif
        enddo
        !**************************************!
        !*** end of iterations             ****!
        !**************************************!

        ! compute final value of new state vector at new time level
        ! result will be in hvy_block again (inout)
        u = uold
        do iter = 1, M_iter+1
            ! call checknan_real( rhs(:,:,:,2,iter), "some rhs " )
            u = u + beta * rhs(:,:,:,:,iter) * phiMat(iter,M_iter+1)
        enddo

        call fft3(inx=u, outk=uk)

        call checknan_real( u(:,:,:,2), "did krylov " )

        if (mpirank==0) then
            open(14,file='krylov_err.t',status='unknown',position='append')
            write (14,'(3(g15.8,1x),2(i3,1x))') time, dt, err, M_iter, M_max
            close(14)
        endif

    end subroutine krylov


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function expM_pade(H) result(E)
        implicit none
        real(kind=pr),dimension(:,:),intent(in)                 :: H
        real(kind=pr),dimension(size(H,1),size(H,2))            :: E

        real(kind=pr),dimension(size(H,1),size(H,2))            :: expH
        integer                                                 :: m,ideg,lwsp,iexph,ldh,i,j,ii,iflag,ns
        real(kind=pr),allocatable,dimension(:)                  :: wsp
        integer,allocatable,dimension(:)                        :: ipiv
        integer::k

        !*** Berechnung des Matrixexponenten **!
        m     = size(H,1)
        ideg  = 6
        lwsp  = 4*m*m+ideg+1
        iexph = 1
        ldh   = m
        allocate(ipiv(m))   ! integer
        allocate(wsp(lwsp)) ! rk

        !                  !1.0_pr!, da dt*H -> H uebergeben wird, sonst steht hier dt
        call  DGPADM(ideg,m,1.0_pr,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )
        if(iflag.lt.0) stop "error in computing exp(t*H)"

        do j=1,m
            do i=1,m
                ii = (i+iexph-1) + (j-1)*m
                E(i,j) = wsp(ii)
            end do
        end do
        !**************************************!

    end function expM_pade


    !  *----------------------------------------------------------------------|
    subroutine DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )
        implicit none
        integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
        double precision t, H(ldh,m), wsp(lwsp)

        !    *-----Purpose----------------------------------------------------------|
        !    *
        !    *     Computes exp(t*H), the matrix exponential of a general matrix in
        !    *     full, using the irreducible rational Pade approximation to the
        !    *     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
        !      *     combined with scaling-and-squaring.
        !      *
        !      *-----Arguments--------------------------------------------------------|
        !      *
        !      *     ideg      : (input) the degre of the diagonal Pade to be used.
        !      *                 a value of 6 is generally satisfactory.
        !      *
        !      *     m         : (input) order of H.
        !      *
        !*     H(ldh,m)  : (input) argument matrix.
        !*
        !*     t         : (input) time-scale (can be < 0).
        !*
        !*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
        !*
        !*     ipiv(m)   : (workspace)
        !*
        !*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
        !*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
        !*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*                 NOTE: if the routine was called with wsp(iptr),
        !*                       then exp(tH) will start at wsp(iptr+iexph-1).
        !*
        !*     ns        : (output) number of scaling-squaring used.
        !*
        !*     iflag     : (output) exit flag.
        !*                      0 - no problem
        !*                     <0 - problem
        !*
        !*----------------------------------------------------------------------|
        !*     Roger B. Sidje (rbs@maths.uq.edu.au)
        !*     EXPOKIT: Software Package for Computing Matrix Exponentials.
        !*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
        !*----------------------------------------------------------------------|
        !*
        integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
        double precision hnorm,scale,scale2,cp,cq

        intrinsic INT,ABS,DBLE,LOG,MAX

        !*---  check restrictions on input parameters ...
        mm = m*m
        iflag = 0
        if ( ldh.lt.m ) iflag = -1
        if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
        if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPADM)'
        !*
        !*---  initialise pointers ...
        !*
        icoef = 1
        ih2 = icoef + (ideg+1)
        ip  = ih2 + mm
        iq  = ip + mm
        ifree = iq + mm
        !*
        !*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
        !*     and set scale = t/2^ns ...
        !*
        do i = 1,m
            wsp(i) = 0.0d0
        enddo
        do j = 1,m
            do i = 1,m
                wsp(i) = wsp(i) + ABS( H(i,j) )
            enddo
        enddo
        hnorm = 0.0d0
        do i = 1,m
            hnorm = MAX( hnorm,wsp(i) )
        enddo
        hnorm = ABS( t*hnorm )

        ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
        scale = t / DBLE(2**ns)
        scale2 = scale*scale
        !*
        !*---  compute Pade coefficients ...
        !*
        i = ideg+1
        j = 2*ideg+1
        wsp(icoef) = 1.0d0
        do k = 1,ideg
            wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
        enddo
        !*
        !*---  H2 = scale2*H*H ...
        !*
        call DGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
        !*
        !*---  initialize p (numerator) and q (denominator) ...
        !*
        cp = wsp(icoef+ideg-1)
        cq = wsp(icoef+ideg)
        do j = 1,m
            do i = 1,m
                wsp(ip + (j-1)*m + i-1) = 0.0d0
                wsp(iq + (j-1)*m + i-1) = 0.0d0
            enddo
            wsp(ip + (j-1)*(m+1)) = cp
            wsp(iq + (j-1)*(m+1)) = cq
        enddo
        !*
        !*---  Apply Horner rule ...
        !*
        iodd = 1
        k = ideg - 1
        100  continue
        iused = iodd*iq + (1-iodd)*ip
        call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iused),m, &
        wsp(ih2),m, 0.0d0,wsp(ifree),m )
        do j = 1,m
            wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
        enddo
        ip = (1-iodd)*ifree + iodd*ip
        iq = iodd*ifree + (1-iodd)*iq
        ifree = iused
        iodd = 1-iodd
        k = k-1
        if ( k.gt.0 )  goto 100
        !*
        !*---  Obtain (+/-)(I + 2*(p\q)) ...
        !*
        if ( iodd .eq. 1 ) then
            call DGEMM( 'n','n',m,m,m, scale,wsp(iq),m, &
            H,ldh, 0.0d0,wsp(ifree),m )
            iq = ifree
        else
            call DGEMM( 'n','n',m,m,m, scale,wsp(ip),m, &
            H,ldh, 0.0d0,wsp(ifree),m )
            ip = ifree
        endif
        call DAXPY( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
        call DGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
        if ( iflag.ne.0 ) stop 'Problem in DGESV (within DGPADM)'
        call DSCAL( mm, 2.0d0, wsp(ip), 1 )
        do j = 1,m
            wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
        enddo
        iput = ip
        if ( ns.eq.0 .and. iodd.eq.1 ) then
            call DSCAL( mm, -1.0d0, wsp(ip), 1 )
            goto 200
        endif
        !*
        !*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
        !*
        iodd = 1
        do k = 1,ns
            iget = iodd*ip + (1-iodd)*iq
            iput = (1-iodd)*ip + iodd*iq
            call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m, &
            0.0d0,wsp(iput),m )
            iodd = 1-iodd
        enddo
        200  continue
        iexph = iput
    END subroutine DGPADM
    !*----------------------------------------------------------------------|



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_sum_all(inout)
        use vars
        implicit none
        real(kind=pr),intent(inout)  :: inout
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(kind=pr) :: tmp
        integer :: mpierr

        call MPI_ALLREDUCE(inout,tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
        inout = tmp

    end subroutine get_sum_all


    subroutine kry_norm(field, norm)
        use vars
        implicit none
        !> heavy data array - block data
        real(kind=pr), intent(inout) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
        !> this is the output of the function:
        real(kind=pr), intent(out)   :: norm

        integer :: ix, iy, iz

        norm = 0.0_pr

        do iz=ra(3),rb(3)
            do iy=ra(2),rb(2)
                do ix=ra(1),rb(1)
                    norm = norm + sum( field(ix,iy,iz,:)**2 )
                enddo
            enddo
        enddo
        call get_sum_all(norm)
        norm = sqrt(norm)

    end subroutine kry_norm

    subroutine scalarproduct(field1, field2, result)
        use vars
        implicit none
        real(kind=pr), intent(inout) :: field1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
        real(kind=pr), intent(inout) :: field2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
        real(kind=pr), intent(out)   :: result

        integer :: ix, iy, iz, ieqn

        result = 0.0_pr

        do iz=ra(3),rb(3)
            do iy=ra(2),rb(2)
                do ix=ra(1),rb(1)
                    do ieqn = 1, neq
                        result = result + field1(ix,iy,iz,ieqn)*field2(ix,iy,iz,ieqn)
                    enddo
                enddo
            enddo
        enddo

        call get_sum_all(result)
    end subroutine scalarproduct

end module krylov_module
