!*********************************************************************!
!        calculate the absorption spectrum oscillator strength       !
!*********************************************************************!
subroutine abs_osc()
    use common_variables
    implicit none
#include<petsc/finclude/petscdef.h>
    PetscInt:: stateE, hx
    PetscScalar, pointer :: gr_v(:)
    PetscScalar, pointer :: wr_v(:)
    PetscInt:: Pstart, Pend
    PetscInt:: stateG
    real*8:: mab_osc_x(nev)
    real*8:: mab_osc_y(nev)
    real*8:: mab_osc_z(nev)
    
    !create the eigenvector with proper dimensions
    call MatCreateVecs( h_sparse, wr, wi, ierr )
    call MatCreateVecs( h_sparse, gr, gi, ierr )
    call VecGetOwnershipRange(gr , Pstart, Pend, ierr)

    if ( rank==0 ) print*, 'Calculating the absorption oscillator strength'
    if(.not. allocated(ab_osc_x)) allocate (ab_osc_x(nev))
    if(.not. allocated(ab_osc_y)) allocate (ab_osc_y(nev))
    if(.not. allocated(ab_osc_z)) allocate (ab_osc_z(nev))
    if(.not. allocated(ab_sys_eval)) allocate (ab_sys_eval(nev))
    ab_osc_x = 0.d0
    ab_osc_y = 0.d0
    ab_osc_z = 0.d0
    ab_sys_eval = 0.d0
    
    !select the ground state
    stateG = 1
    call EPSGetEigenpair(eps_cxt_sys, stateG-1, ger, gei, gr, gi, ierr )
    call VecGetArrayReadF90(gr,gr_v,ierr)
    
    !loop over excited states
    do stateE = 2, nev
        call EPSGetEigenpair(eps_cxt_sys, stateE-1, er, ei, wr, wi, ierr )
        call VecGetArrayReadF90(wr,wr_v,ierr)

        !transition energy
        ab_sys_eval(stateE) = er - ger
        ! print*, ab_sys_eval(stateE)
        do hx = Pstart+1, Pend
            !the oscillator strength
            ab_osc_x(stateE) = ab_osc_x(stateE) + gr_v(hx - Pstart) * ux(hx) * wr_v(hx - Pstart)
            ab_osc_y(stateE) = ab_osc_y(stateE) + gr_v(hx - Pstart) * uy(hx) * wr_v(hx - Pstart)
            ab_osc_z(stateE) = ab_osc_z(stateE) + gr_v(hx - Pstart) * uz(hx) * wr_v(hx - Pstart)
        end do
        
        call VecRestoreArrayReadF90(wr,wr_v,ierr)
    end do
    call VecRestoreArrayReadF90(gr,gr_v,ierr)

!   print*, Pstart
!   call MPI_REDUCE(Pstart, sum_Pstart, 1, MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
!   print*, sum_Pstart

    call VecDestroy( wr, ierr )
    call VecDestroy( wi, ierr )
    call VecDestroy( gr, ierr )
    call VecDestroy( gi, ierr )
    
    call MPI_Reduce(ab_osc_x, mab_osc_x, nev, MPI_REAL8, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
!    call MPI_Bcast(mab_osc_x, nev, MPI_REAL8,0,MPI_COMM_WORLD, ierr)
    ab_osc_x =  mab_osc_x
    call MPI_Reduce(ab_osc_y, mab_osc_y, nev, MPI_REAL8, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    ab_osc_y = mab_osc_y
    call MPI_Reduce(ab_osc_z, mab_osc_z, nev, MPI_REAL8, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    ab_osc_z = mab_osc_z
    
    if(rank == 0) then
        ab_osc_x = ab_osc_x ** 2
        ab_osc_y = ab_osc_y ** 2
        ab_osc_z = ab_osc_z ** 2
    end if
    
    ! print*, ab_osc_x
    ! print*, ab_osc_y
    ! print*, ab_osc_z
    ! do stateE = 1, nev
        ! print*, ab_osc_x(stateE), ab_osc_y(stateE),ab_osc_z(stateE)
    ! end do
end subroutine