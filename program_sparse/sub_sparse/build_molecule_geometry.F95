!*******************************************************************************!
!   build the geometry of molecule
!*******************************************************************************!
subroutine build_molecule_geometry()
    use common_variables
    implicit none
    
    integer n, da
    ! integer n1, n2
    !x: long axis 
    !y: short axis
    !z: pi stack

    !dipolar with three dimensions
    ! if (.not. allocated(molecule_pos)) allocate(molecule_pos(nmax,es_zmax,3))
    ! molecule_pos = 0.d0
    
!    do n = 1, nmax
!    do da = 1, 3 !quadrupolar
!        molecule_pos(n,da,1) = (n-1)*x_offset + (da -1.0d0) * armlength
!        molecule_pos(n,da,2) = (n-1)*y_offset
!        molecule_pos(n,da,3) = (n-1)*z_offset
!    end do
!    end do
    
!    print*, dcos(pi/2.d0)
!    do n = 1, nmax
!    do da = 1, 3 !quadrupolar
!            print*, dcos((n-1)*agg_angle) 
!        molecule_pos(n,da,1) = (-1)**(da-1)*armlength/2.d0*dcos((n-1)*agg_angle) 
!        molecule_pos(n,da,2) = (-1)**(da-1)*armlength/2.d0*dsin((n-1)*agg_angle) 
!        molecule_pos(n,da,3) = (n-1)*z_offset
!    end do
!    end do

    ! ref_pos()

    if (.not. periodic) then
        !dipolar with three dimensions
        if (.not. allocated(molecule_pos)) allocate(molecule_pos(nmax,es_zmax,3))
        molecule_pos = 0.d0

        do n = 1, nmax
        do da = 1, 3 !quadrupolar
            ! molecule_pos(n,da,1) = (n-1)*x_offset + (da -1.0d0) * armlength * (-1)**(n-1)
            molecule_pos(n,da,1) = (n-1)*x_offset + (da -1.0d0) * armlength 
            molecule_pos(n,da,2) = (n-1)*y_offset
            molecule_pos(n,da,3) = (n-1)*z_offset
        end do
        end do
    else if (periodic .and. nearest_neighbor) then
        !dipolar with three dimensions
        if (.not. allocated(molecule_pos)) allocate(molecule_pos(0:nmax+1,es_zmax,3))
        molecule_pos = 0.d0

        do n = 0, nmax + 1
        do da = 1, 3 !quadrupolar
            ! molecule_pos(n,da,1) = (n-1)*x_offset + (da -1.0d0) * armlength * (-1)**(n-1)
            molecule_pos(n,da,1) = (n-1)*x_offset + (da -1.0d0) * armlength 
            molecule_pos(n,da,2) = (n-1)*y_offset
            molecule_pos(n,da,3) = (n-1)*z_offset
        end do
        end do
    else if (periodic .and. (.not. nearest_neighbor)) then
        !dipolar with three dimensions
        if (.not. allocated(molecule_pos)) allocate(molecule_pos(2-nmax/2:nmax*3/2,es_zmax,3))
        molecule_pos = 0.d0

        do n = 2 - nmax/2,nmax*3/2
        do da = 1, 3 !quadrupolar
            ! molecule_pos(n,da,1) = (n-1)*x_offset + (da -1.0d0) * armlength * (-1)**(n-1)
            molecule_pos(n,da,1) = (n-1)*x_offset + (da - 1.0d0) * armlength 
            molecule_pos(n,da,2) = (n-1)*y_offset
            molecule_pos(n,da,3) = (n-1)*z_offset
            ! print*, molecule_pos(n,da,1), molecule_pos(n,da,2), molecule_pos(n,da,3)
            ! print*, molecule_pos(1,2,1)
        end do
        end do
    end if

        !allocate the charge
    if(.not. allocated(tq)) allocate(tq(es_zmax, 3))
    
    tq = 0.d0
    tq(es_z1,leftdonor) = 1.d0
    tq(es_z1,acceptor) = -1.d0
    tq(es_z2,rightdonor) = 1.d0
    tq(es_z2,acceptor) = -1.d0
    
end subroutine