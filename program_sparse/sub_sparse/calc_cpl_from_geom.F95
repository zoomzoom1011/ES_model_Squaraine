!*********************************************************************!
!    Calculate the coulomb couplings from geometry
!*********************************************************************!
subroutine calc_cpl_from_geom()
    use common_variables
    implicit none

#include<petsc/finclude/petscdef.h>
    integer, allocatable :: atype(:,:)
    integer tmp2, anum2, i, j, &
            sub, la, lb, lc, labc, sub2, da, db, dc, dabc, &
            tqf_no, fno, anum
    character*100 read_line, tmp
    integer:: da1, da2, dn, n1, n2, n3, n
    integer:: state1, state2
    PetscInt:: sysnx
    !constants for au to cm-1 conversion
    real*8, parameter :: eo         = 8.854187817d-12 !(f/m)
    real*8, parameter :: plancks    = 6.62606957d-34  !kg m**2 s**-2
    real*8, parameter :: csol       = 2.99792458d8    !m s**-1
    real*8, parameter :: aucharge_c = 1.602176565d-19 !c per au
    real*8, parameter :: au_debye   = 2.54175d0       !au per debye
    real*8:: coeff
    real*8:: dmat(3)
    integer:: get_charge
    real*8:: get_distance
    ! real*8, allocatable:: mol2pos(:, :, :)
    integer:: bring_inside_aggregate_range

    ! allocate(mol2pos(kount_duplicate_stacks, 2, 3))
    
    anum = 3
    !=================================================================!
    !    calculate the coulomb coupling based on the tdm and geometry
    !=================================================================!
    if (.not. allocated(coulomb_coupling)) allocate(coulomb_coupling(sys_kount))
    coulomb_coupling = 0.d0
    
    if(nmax > 1) then
        do sysnx = 1, sys_kount
            if (periodic .and. nearest_neighbor) then 
                do n1 = 1, nmax
                    state1 = sys_state(sysnx,n1)%es_state
                    !periodic and nearest neighbor
                    do dn = -1, 1
                        if (dn == 0) cycle 
                        n2 = n1 + dn
                        n3 = bring_inside_aggregate_range(n2)
                        state2 = sys_state(sysnx,n3)%es_state
                        
                        do da1 = 1, anum
                        do da2 = 1, anum
                            coulomb_coupling(sysnx) =coulomb_coupling(sysnx)+ &
                                     get_charge(state1, da1) &
                                    *get_charge(state2, da2) &
                                    /get_distance(n1,da1,n2,da2)
                        end do
                        end do
                    end do
                end do
                coulomb_coupling(sysnx) = coulomb_coupling(sysnx)/2.d0
                If ( nmax == 2 ) then
                    coulomb_coupling(sysnx) = coulomb_coupling(sysnx)/2.D0        !There should only be one coupling between 1 and 2, above i overcount by 1
                End If
            else if (.not. periodic .and. nearest_neighbor) then 
                do n1 = 1, nmax-1
                    state1 = sys_state(sysnx,n1)%es_state
                    !open boundary and nearest neighbor
                    n2 = n1 + 1
                    state2 = sys_state(sysnx,n2)%es_state
                    
!                    if (state1 == es_z .and. state2 == es_z ) then
!                            coulomb_coupling(sysnx) =coulomb_coupling(sysnx)+ v 
!                    end if 
                    do da1 = 1, anum
                    do da2 = 1, anum
                            coulomb_coupling(sysnx) =coulomb_coupling(sysnx)+ &
                                     get_charge(state1, da1) &
                                    *get_charge(state2, da2) &
                                    /get_distance(n1,da1,n2,da2)
                    end do
                    end do
                end do
            else if ((.not. periodic) .and. (.not. nearest_neighbor)) then 
                do n1 = 1, nmax-1
                    state1 = sys_state(sysnx,n1)%es_state
                    !open boundary and not nearest neighbor
                do n2 = n1+1, nmax
                    state2 = sys_state(sysnx,n2)%es_state
                    do da1 = 1, anum
                    do da2 = 1, anum
                            coulomb_coupling(sysnx) =coulomb_coupling(sysnx)+ &
                                     get_charge(state1, da1) &
                                    *get_charge(state2, da2) &
                                    /get_distance(n1,da1,n2,da2)
                    end do
                    end do
                end do
                end do
            else if(periodic .and. (.not. nearest_neighbor)) then
                do n1 = 1, nmax
                    state1 = sys_state(sysnx,n1)%es_state
                    !periodic and not nearest neighbor
                do dn = -nmax/2 + 1, nmax/2
                    if (dn == 0) cycle 
                    n2 = n1 + dn
                    n3 = bring_inside_aggregate_range(n2)
                    state2 = sys_state(sysnx,n3)%es_state
                    do da1 = 1, anum
                    do da2 = 1, anum
                            coulomb_coupling(sysnx) =coulomb_coupling(sysnx)+ &
                                     get_charge(state1, da1) &
                                    *get_charge(state2, da2) &
                                    /get_distance(n1,da1,n2,da2)
                    end do
                    end do
                end do
                end do
                coulomb_coupling(sysnx) = coulomb_coupling(sysnx)/2.d0
                ! If ( nmax == 2 ) then
                    ! coulomb_coupling(sysnx) = coulomb_coupling(sysnx)/2.D0        !There should only be one coupling between 1 and 2, above i overcount by 1
                ! End If
            end if
            
            !use other stacks to correct the comlomb coupling
            ! print*,'before:',coulomb_coupling(sysnx)
            ! ! if(duplicate_stacks) then

                ! ! do n1 = 1, nmax
                    ! ! state1 = sys_state(sysnx,n1)%es_state
                    ! ! if(state1 == es_n) cycle
                    ! ! n = (n1-1)/2
                    ! ! n2 = 4*n+3-n1

                    ! ! do la = -duplicate_stacks_a, duplicate_stacks_a
                    ! ! do lb = -duplicate_stacks_b, duplicate_stacks_b
                    ! ! do lc = -duplicate_stacks_c, duplicate_stacks_c
                        ! ! if(lb == 0 .and. la == 0) cycle
                        ! ! labc = nx_lattice(la, lb, lc)
                        ! ! ! print*, la, lb, lc, labc
                        ! ! if(labc == empty) cycle
                        ! ! dmat = (/la,lb,lc/)
                        ! ! ! n1
                        ! ! mol2pos(labc,:,:) = mol1pos(n1,:,:)
                        ! ! mol2pos(labc,:,1) = mol2pos(labc,:,1) + dot_product( to_cartesian(1,:), dmat )
                        ! ! mol2pos(labc,:,2) = mol2pos(labc,:,2) + dot_product( to_cartesian(2,:), dmat )
                        ! ! mol2pos(labc,:,3) = mol2pos(labc,:,3) + dot_product( to_cartesian(3,:), dmat )
                        
                        ! ! do da1 = 1, anum
                        ! ! do da2 = 1, anum
                            ! ! distance = (mol2pos(labc, da1, 1)- mol1pos(n1, da2, 1))**2 + &
                                    ! ! (mol2pos(labc, da1, 2)- mol1pos(n1, da2, 2))**2 + &
                                    ! ! (mol2pos(labc, da1, 3)- mol1pos(n1, da2, 3))**2
                            ! ! distance = dsqrt(distance)
                            ! ! coulomb_coupling(sysnx) =coulomb_coupling(sysnx)+ &
                                ! ! (get_charge(es_z, da1)* get_charge(es_z, da2))*1.d0/distance
                        ! ! end do
                        ! ! end do
                        
                        ! ! ! n2
                        ! ! mol2pos(labc,:,:) = mol1pos(n2,:,:)
                        ! ! mol2pos(labc,:,1) = mol2pos(labc,:,1) + dot_product( to_cartesian(1,:), dmat )
                        ! ! mol2pos(labc,:,2) = mol2pos(labc,:,2) + dot_product( to_cartesian(2,:), dmat )
                        ! ! mol2pos(labc,:,3) = mol2pos(labc,:,3) + dot_product( to_cartesian(3,:), dmat )
                        ! ! do da1 = 1, anum
                        ! ! do da2 = 1, anum
                            ! ! distance = (mol2pos(labc, da1, 1)- mol1pos(n1, da2, 1))**2 + &
                                    ! ! (mol2pos(labc, da1, 2)- mol1pos(n1, da2, 2))**2 + &
                                    ! ! (mol2pos(labc, da1, 3)- mol1pos(n1, da2, 3))**2
                            ! ! distance = dsqrt(distance)
                            ! ! coulomb_coupling(sysnx) =coulomb_coupling(sysnx)+ &
                                ! ! (get_charge(es_z, da1)* get_charge(es_z, da2))*1.d0/distance
                        ! ! end do
                        ! ! end do
                    ! ! end do
                    ! ! end do
                    ! ! end do
                ! ! end do
            ! ! end if
            ! print*,'after:',coulomb_coupling(sysnx)
            ! n1 = 0
            ! do la = 1, UC
            ! do sub = 1, MPUC
                ! n1 = n1 + 1
                ! state1 = sys_state(sysnx,n1)%es_state
                ! refpos(:,:) = mol1pos(n1,:,:)
                
                ! n2 = 0
                ! do la1 = 1, UC
                ! do sub1 = 1, MPUC
                    ! n2 = n2 + 1
                    ! if(n2 <= n1) cycle
                    ! state2 = sys_state(sysnx,n2)%es_state
                    
                ! do lb1 = -duplicate_stacks_b, duplicate_stacks_b
                ! do lc1 = -duplicate_stacks_c, duplicate_stacks_c
                    ! if(la1 == la .and. sub1 == sub .and. lb1 == 0 .and. lc1 == 0) cycle
                    ! dmat = (/la1-la,lb1,lc1/)
                    
                    !!calculate the coordinates of the second molecule
                    ! if ( .not. allocated( mol2pos ) ) allocate( mol2pos( anum, 3 ) )
                    ! mol2pos(:,:) = refpos_cartesian( sub1, :, : )
                    ! mol2pos(:,1) = mol2pos(:,1) + dot_product( to_cartesian(1,:), dmat )
                    ! mol2pos(:,2) = mol2pos(:,2) + dot_product( to_cartesian(2,:), dmat )
                    ! mol2pos(:,3) = mol2pos(:,3) + dot_product( to_cartesian(3,:), dmat )
                    
                    !!go over all tq in each molecule and calculate the coupling
                    ! do anum1 = 1, anum
                    ! do anum2 = 1, anum
                        ! coulomb_coupling(sysnx) =                            &
                            ! coulomb_coupling(sysnx) +                        &
                            ! get_charge(state1,anum1)*get_charge(sub2,anum2)/dsqrt(              &
                            ! ( mol2pos(anum1,1) - refpos(anum2,1) ) ** 2                    &
                           ! +( mol2pos(anum1,2) - refpos(anum2,2) ) ** 2                    &
                           ! +( mol2pos(anum1,3) - refpos(anum2,3) ) ** 2 )
                    ! end do
                    ! end do
                ! end do
                ! end do
                    
                ! end do
                ! end do
            ! end do
            ! end do
        end do
    end if

    coeff = (aucharge_c**2)*(1.d10)/(4.d0*pi*100.d0*Eo*plancks*csol)/dielectric
    coulomb_coupling = coeff * coulomb_coupling

    if (no_frenkel) coulomb_coupling = 0.d0
    ! if ( rank==0 ) print*, 'max coulomb coupling:', coulomb_coupling(sys_kount)
!    print*, coulomb_coupling*hw
    
    if (.not. allocated(ux) ) allocate(ux(sys_kount))
    if (.not. allocated(uy) ) allocate(uy(sys_kount))
    if (.not. allocated(uz) ) allocate(uz(sys_kount))
    ux = 0.d0
    uy = 0.d0
    uz = 0.d0
    
    do sysnx = 1, sys_kount
        do n1 = 1, nmax
            state1 = sys_state(sysnx,n1)%es_state
            do da1 = 1, anum
                ux(sysnx) = ux(sysnx) + get_charge(state1, da1)* mol1pos(n1, da1, 1)
                uy(sysnx) = uy(sysnx) + get_charge(state1, da1)* mol1pos(n1, da1, 2)
                uz(sysnx) = uz(sysnx) + get_charge(state1, da1)* mol1pos(n1, da1, 3)
            end do
        end do 
!         print*, ux(sysnx),uy(sysnx),uz(sysnx)
    end do
    
end subroutine
!*****************************************************************************!
real*8 function get_distance(n1,da1,n2,da2)
    use common_variables
    implicit none

    integer, intent(in) :: n1, n2, da1, da2

    get_distance = (mol1pos(n1, da1, 1)-mol1pos(n2, da2, 1))**2&
                  +(mol1pos(n1, da1, 2)-mol1pos(n2, da2, 2))**2&
                  +(mol1pos(n1, da1, 3)-mol1pos(n2, da2, 3))**2
    get_distance = dsqrt(get_distance)
    ! !convert to meters from angstroms
    ! get_distance = get_distance * 1.d-10
    return
end function
!*****************************************************************************!
integer function get_charge( state, da )
    use common_variables
    implicit none

    integer, intent(in) :: state, da

    get_charge = 0
    if ( state == es_n ) then
        continue
    else if ( state == es_z1 ) then
        if ( da == leftdonor ) then
            get_charge = 1
        else if ( da == acceptor ) then
            get_charge = -1
        end if
    else if ( state == es_z2 ) then
        if ( da == rightdonor ) then
            get_charge = 1
        else if ( da == acceptor ) then
            get_charge = -1
        end if
    end if
    return
end function
