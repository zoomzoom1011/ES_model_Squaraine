!*********************************************************************!
!        calculate coulomb coupling of the dimer zwitter state  !
!*********************************************************************!
subroutine  calc_coupling()
    use common_variables
    implicit none
        
    integer n, da

    real*8:: distance
    integer:: n1, n2, da1, da2, sysnx
    integer:: state1, state2
    real*8:: coeff
    
    !physical constants
    !constants for au to cm-1 conversion
    real*8, parameter :: eo = 8.854187817d-12 !(f/m)
    real*8, parameter :: plancks  = 6.62606957d-34  !kg m**2 s**-2
    real*8, parameter :: csol  = 2.99792458d8    !m s**-1
    real*8, parameter :: aucharge_c  = 1.602176565d-19 !c per au
    
    !constants for au to cm-1 conversion
!    real*8, parameter :: EO = 8.854187817620d-12 !(f/m)
!    real*8, parameter :: PLANCKS  = 6.62606957d-34  !kg m**2 s**-2
!    real*8, parameter :: CSOL  = 2.99792458d8    !m s**-1
!    real*8, parameter :: aucharge_c  = 1.602176565d-19 !c per au
!    real*8, parameter :: au_debye      = 2.54175d0
!    real*8, parameter :: joule_ev      = 6.24150934d18
!    real*8 tq(3,3)
!    integer:: anum = 3
!    integer m, l
!    integer i, j
!    integer anum2
!    integer lx, ly, lz, lxyz
    
    !allocate the charge according to the model
    !here, the armlength is 5A, offset in z axis is 3.5A
    !left is cation, right is left
    !first molecule locate in the center of cartesian coordinate
!    if (.not. allocated(molpos_cartesian)) allocate(molpos_cartesian(kount_lattice,anum,3))
!    molpos_cartesian = 0.d0

    !dipolar with three dimensions
    if (.not. allocated(molecule_pos)) allocate(molecule_pos(nmax,es_zmax,3))
    molecule_pos = 0.d0

    do n = 1, nmax
    do da = 1, 3 !quadrupole
        molecule_pos(n,da,1) = (n-1)*x_offset + (da -2.0d0) * armlength
        molecule_pos(n,da,2) = (n-1)*y_offset
        molecule_pos(n,da,3) = (n-1)*z_offset
    end do
    end do

    !allocate the charge
    if(.not. allocated(tq)) allocate(tq(es_zmax, 3))
    
    tq = 0.d0
    tq(es_z1,leftdonor) = 1.d0
    tq(es_z1,acceptor) = -1.d0
    tq(es_z2,rightdonor) = 1.d0
    tq(es_z2,acceptor) = -1.d0
    
    if (.not. allocated(coulomb_coupling)) allocate(coulomb_coupling(sys_kount))
    coulomb_coupling = 0.d0
    
    if (nmax > 1) then
        do sysnx = 1, sys_kount
            do n1 = 1, nmax-1
                state1 = sys_state(sysnx,n1)%es_state
            do n2 = n1+1, nmax
                state2 = sys_state(sysnx,n2)%es_state
                
                do da1 = 1, 3
                do da2 = 1, 3
                    distance = (molecule_pos(n1, da1, 1)- molecule_pos(n2, da2, 1))**2 + &
                                (molecule_pos(n1, da1, 2)- molecule_pos(n2, da2, 2))**2 + &
                                (molecule_pos(n1, da1, 3)- molecule_pos(n2, da2, 3))**2
                    distance = dsqrt(distance)
                    coulomb_coupling(sysnx) = coulomb_coupling(sysnx)+ tq(state1, da1)* tq(state2, da2)/distance
                        
                end do
                end do
            end do
            end do
!            print*, coulomb_coupling(sysnx)
        end do
    end if
    
    
    coeff = (aucharge_c**2)*(1.d10)/(4.d0*pi*100.d0*Eo*plancks*csol)/hw/dielectric

    coulomb_coupling = coeff * coulomb_coupling
    !scalar
    coulomb_coupling = coulomb_coupling/(nmax-1)

    if (no_frenkel) coulomb_coupling = 0.d0
    
    ! print*, coulomb_coupling*hw
    
    if (.not. allocated(ux) ) allocate(ux(sys_kount))
    if (.not. allocated(uy) ) allocate(uy(sys_kount))
    if (.not. allocated(uz) ) allocate(uz(sys_kount))
    ux = 0.d0
    uy = 0.d0
    uz = 0.d0
    
    do sysnx = 1, sys_kount
        do n1 = 1, nmax
            state1 = sys_state(sysnx,n1)%es_state
            do da1 = 1, 3
                ux(sysnx) = ux(sysnx) + tq(state1, da1)* molecule_pos(n1, da1, 1)
                uy(sysnx) = uy(sysnx) + tq(state1, da1)* molecule_pos(n1, da1, 2)
                uz(sysnx) = uz(sysnx) + tq(state1, da1)* molecule_pos(n1, da1, 3)
            end do
        end do
    end do
    
    
  
!    molpos_cartesian(1, 1, 1) = 0.d0
!    molpos_cartesian(1, 1, 2) = 0.d0
!    molpos_cartesian(1, 1, 3) = 0.d0
!    molpos_cartesian(1, 2, 1) = armlength
!    molpos_cartesian(1, 2, 2) = 0.d0
!    molpos_cartesian(1, 2, 3) = 0.d0
!    molpos_cartesian(1, 3, 1) = -armlength
!    molpos_cartesian(1, 3, 2) = 0.d0
!    molpos_cartesian(1, 3, 3) = 0.d0
!
!    tq(1) = 1.d0
!    tq(2) = -1.d0
!    do lx = 1, lattice_x
!    do ly = 1, lattice_y
!    do lz = 1, lattice_z
!        lxyz = nx_lattice(lx, ly, lz)
!    do anum2 = 1, anum
!        molpos_cartesian( lxyz, anum2, : ) =        &
!            molpos_cartesian( 1, anum2, : )+        &
!            (/(lx-1)*x_offset,(ly-1)*y_offset,(lz-1)*z_offset/)
!            print*, molpos_cartesian( lxyz, anum2, : )
!    end do
!    end do
!    end do
!    end do
!
!    
!    
!    if (.not. allocated(coulomb_coupling)) allocate(coulomb_coupling(nmax, nmax))
!    coulomb_coupling = 0.d0
!    
!    do l = 1, nmax
!    do m = 1, nmax
        !no self-coupling
!        if ( l >= m) cycle
        !go over all tq
!        do i = 1, anum
!        do j = 1, anum
            !j = tq1*tq2/R
!            coulomb_coupling(l, m) = coulomb_coupling(l, m) +   &
!                         tq(i)*tq(j)/dsqrt(                          &
!                   (molpos_cartesian(l,i,1)-molpos_cartesian(m,j,1))**2    &  
!               + (molpos_cartesian(l,i,2)-molpos_cartesian(m,j,2))**2    & 
!               + (molpos_cartesian(l,i,3)-molpos_cartesian(m,j,3))**2 )
!        end do
!        end do
!    end do
!    end do
!    
!    coulomb_coupling = coulomb_coupling *  (aucharge_c**2)*(1.d10)/(4.d0*pi*Eo)/hw*joule_ev*ev
!    print*, coulomb_coupling
!        
!        if (.not. allocated(molecule_posx)) allocate(molecule_posx(kount_lattice))
!    if (.not. allocated(molecule_posy)) allocate(molecule_posy(kount_lattice))
!    if (.not. allocated(molecule_posz)) allocate(molecule_posz(kount_lattice))
!    molecule_posx = 0.d0
!    molecule_posy = 0.d0
!    molecule_posz = 0.d0
!    
!    do lx = 1, lattice_x
!    do ly = 1, lattice_y
!    do lz = 1, lattice_z
!        lxyz = nx_lattice(lx, ly, lz)
!        molecule_posx(lxyz) = sum(molpos_cartesian(lxyz,:,1))/(1.d0*anum)
!        molecule_posy(lxyz) = sum(molpos_cartesian(lxyz,:,2))/(1.d0*anum)
!        molecule_posz(lxyz) = sum(molpos_cartesian(lxyz,:,3))/(1.d0*anum)
!    end do
!    end do
!    end do
end subroutine