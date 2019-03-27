!*********************************************************************!
!                pl spectrum to each state         !
!*********************************************************************!
subroutine pl_osc()
    use common_variables
    implicit none
    
#include<petsc/finclude/petscdef.h>
    PetscInt:: statelo, stateup, hx
    PetscScalar, pointer :: gr_v(:)
    PetscScalar, pointer :: wr_v(:)
    real*8:: mpl_osc_x(nev,nev)
    real*8:: mpl_osc_y(nev,nev)
    real*8:: mpl_osc_z(nev,nev)
    PetscInt:: Pstart, Pend
    
!    if(.not. allocated(sys_pl_osc)) allocate(sys_pl_osc(nev, nev))
!    sys_pl_osc = 0.d0
    
!    do statelo = 1, nev-1
!    do stateup = statelo+1, nev
!        do state = 1, nev
!            sys_pl_osc(statelo,stateup) = sys_pl_osc(statelo,stateup)+ sys_h(state, statelo)*sys_h(state, stateup)*sys_mu(state,state)
!        end do
!            sys_pl_osc(statelo,stateup) = sys_pl_osc(statelo,stateup) * sys_pl_osc(statelo,stateup) 
!    end do
!    end do

    !create the eigenvector with proper dimensions
    call MatCreateVecs( h_sparse, wr, wi, ierr )
    call MatCreateVecs( h_sparse, gr, gi, ierr )
    call VecGetOwnershipRange(gr , Pstart, Pend, ierr)
    
    if ( rank==0 ) print*, 'Calculating the photoluminescence oscillator strength'
    if(.not. allocated(pl_osc_x)) allocate (pl_osc_x(nev, nev))
    if(.not. allocated(pl_osc_y)) allocate (pl_osc_y(nev, nev))
    if(.not. allocated(pl_osc_z)) allocate (pl_osc_z(nev, nev))
    if(.not. allocated(pl_sys_eval)) allocate (pl_sys_eval(nev, nev))
    pl_osc_x = 0.d0
    pl_osc_y = 0.d0
    pl_osc_z = 0.d0
    pl_sys_eval = 0.d0
    
    do statelo = 1, nev-1
        call EPSGetEigenpair(eps_cxt_sys, statelo-1, ger, gei, gr, gi, ierr )
        call VecGetArrayReadF90(gr,gr_v,ierr)
    do stateup = statelo + 1, nev
        call EPSGetEigenpair(eps_cxt_sys, stateup-1, er, ei, wr, wi, ierr )
        call VecGetArrayReadF90(wr,wr_v,ierr)
        
        !transition energy
        pl_sys_eval(statelo,stateup) = er - ger
        
        do hx = Pstart+1, Pend
            pl_osc_x(statelo,stateup) = pl_osc_x(statelo,stateup) + gr_v(hx - Pstart) * ux(hx) * wr_v(hx - Pstart)
            pl_osc_y(statelo,stateup) = pl_osc_y(statelo,stateup) + gr_v(hx - Pstart) * uy(hx) * wr_v(hx - Pstart)
            pl_osc_z(statelo,stateup) = pl_osc_z(statelo,stateup) + gr_v(hx - Pstart) * uz(hx) * wr_v(hx - Pstart)
        end do
        call VecRestoreArrayReadF90(wr,wr_v,ierr)
    end do
    call VecRestoreArrayReadF90(gr,gr_v,ierr)
    end do
    
    call VecDestroy( wr, ierr )
    call VecDestroy( wi, ierr )
    call VecDestroy( gr, ierr )
    call VecDestroy( gi, ierr )
    
    call MPI_Reduce(pl_osc_x, mpl_osc_x, nev*nev, MPI_REAL8, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    pl_osc_x =  mpl_osc_x
    call MPI_Reduce(pl_osc_y, mpl_osc_y, nev*nev, MPI_REAL8, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    pl_osc_y = mpl_osc_y
    call MPI_Reduce(pl_osc_z, mpl_osc_z, nev*nev, MPI_REAL8, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    pl_osc_z = mpl_osc_z
    
    if(rank == 0) then
        pl_osc_x = pl_osc_x ** 2
        pl_osc_y = pl_osc_y ** 2
        pl_osc_z = pl_osc_z ** 2
    end if
    
end subroutine