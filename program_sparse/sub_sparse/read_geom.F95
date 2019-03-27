subroutine read_geom()
    use common_variables
    implicit none
    
    integer, allocatable :: atype(:,:)
    integer anum, tmp2, &
            anum2,&
            sub, la, lb, lc, labc, sub2, da, db, dc, dabc, &
            tqf_no, fno
    character*100 read_line, tmp
    integer:: da1, da2, dn, n1, n2
    integer:: state1, state2, sysnx
    integer:: UC
    real*8:: dmat(3), rotax(3), rotangle, &
           norm, rmat(3,3), a1(3), a2(3), a3(3)
    real*8, allocatable :: refpos_cartesian(:,:,:)
    character*128:: cwd
    Integer:: istat, getcwd

    !Get the current working directory
    istat = getcwd( cwd )
    if ( istat .ne. 0 ) stop 'getcwd: error'

    !Adjust the geometry file name
    do sub = 1, sublattice
        xyz_file(sub) = trim(cwd)//trim(xyz_file(sub))
    end do

    anum = 3
    UC = nmax/sublattice

    ! if ( .not. allocated(tq)    ) allocate(tq(es_zmax,anum))
    if ( .not. allocated(atype) ) allocate(atype(sublattice,anum))
    if ( .not. allocated(refpos_cartesian)) allocate(refpos_cartesian(sublattice,anum,3))
    if ( .not. allocated(mol1pos)) allocate(mol1pos(nmax,anum,3))

    ! ! allocate the charge
    ! tq = 0.d0

    !==================================================!
    !set the reference molecules
    tqf_no = 51
    do sub = 1, sublattice
        !read in the tq for sub1
        if ( rank==0 ) print*, 'Reading the xyz file:', trim(xyz_file(sub))
        open( unit = tqf_no, file = trim(xyz_file(sub)), status = 'old', action = 'read' )
        !read the header, then atom number
        do
            read( tqf_no, *, end=102 ) tmp2, anum2,  &
                   refpos_cartesian( sub, anum2, 1 ),& !x
                   refpos_cartesian( sub, anum2, 2 ),& !y
                   refpos_cartesian( sub, anum2, 3 )   !z
            atype(sub,anum2) = tmp2
        end do
    102 close(tqf_no)
    end do
    
!     print*, refpos_cartesian
    ! !allocate the charge
    ! if(es_zmax > 3) then
        ! tq(es_c1,leftdonor) = 1
        ! tq(es_a,acceptor) = -1
    ! end if

    !first rotate dipole moment vector to xy plane
    !rotate the geometry
    !a1 normalize vector
    a1(1) = refpos_cartesian(1,1,1)-refpos_cartesian(1,2,1)
    a1(2) = refpos_cartesian(1,1,2)-refpos_cartesian(1,2,2)
    a1(3) = refpos_cartesian(1,1,3)-refpos_cartesian(1,2,3)
    norm = dsqrt(sum(a1(:)**2))
    a1 = a1/norm

    !here find the origin vector projection in the xy plane and normalize it
    a2(1) = refpos_cartesian(2,1,1)-refpos_cartesian(2,2,1)
    a2(2) = refpos_cartesian(2,1,2)-refpos_cartesian(2,2,2)
    a2(3) = refpos_cartesian(2,1,3)-refpos_cartesian(2,2,3)
!     print*, bproj(1), bproj(2), refpos_cartesian(1,1,3)-refpos_cartesian(1,2,3)
    ! print*, bproj
    norm = dsqrt(sum(a2(:)**2))
    a2 = a2/norm

    !get the rotation axis
    call cross_product( a2, a1, a3 )   
    norm = dsqrt( sum(a3(:)**2) ) 
    a3 = a3/norm

    !get the rotation axis and angle
    call cross_product( a3, (/0.d0,0.d0,1.d0/), rotax )                                                             
    norm = dsqrt( sum(rotax(:)**2) )
    rotax = rotax/norm
    !find the rotation angle
    rotangle = dacos(dot_product(a3,(/0.d0,0.d0,1.d0/)))

    !use the rodrigues' rotation formula to rotate the reference molecules to
    !the lab frame
    call get_rmat( rmat, rotax, rotangle )
    to_cartesian = matmul( rmat, to_cartesian )
    !rotate the reference molecules
    do sub = 1, sublattice
    do anum2 = 1, anum
        refpos_cartesian( sub, anum2, : ) =    &
            matmul( rmat, refpos_cartesian( sub, anum2, : ) )
    end do
    end do

    !second rotate the sub1+sub2 to be x-axis
    !a1 normalize vector
    a1(1) = refpos_cartesian(1,1,1)-refpos_cartesian(1,2,1)
    a1(2) = refpos_cartesian(1,1,2)-refpos_cartesian(1,2,2)
    a1(3) = 0.d0
    norm = dsqrt(sum(a1(:)**2))
    a1 = a1/norm

    !here find the origin vector projection in the xy plane and normalize it
    a2(1) = refpos_cartesian(2,1,1)-refpos_cartesian(2,2,1)
    a2(2) = refpos_cartesian(2,1,2)-refpos_cartesian(2,2,2)
    a2(3) = 0.d0
!     print*, bproj(1), bproj(2), refpos_cartesian(1,1,3)-refpos_cartesian(1,2,3)
    ! print*, bproj
    norm = dsqrt(sum(a2(:)**2))
    a2 = a2/norm

    !get the rotation axis
    a3 = a1 + a2
    norm = dsqrt( sum(a3(:)**2) ) 
    a3 = a3/norm

    !get the rotation axis and angle
    call cross_product( a3, (/1.d0,0.d0,0.d0/), rotax )                                                             
    norm = dsqrt( sum(rotax(:)**2) )
    rotax = rotax/norm
    !find the rotation angle
    rotangle = dacos(dot_product(a3,(/1.d0,0.d0,0.d0/)))

    !use the rodrigues' rotation formula to rotate the reference molecules to
    !the lab frame
    call get_rmat( rmat, rotax, rotangle )
    to_cartesian = matmul( rmat, to_cartesian )
    !rotate the reference molecules
    do sub = 1, sublattice
    do anum2 = 1, anum
        refpos_cartesian( sub, anum2, : ) =    &
            matmul( rmat, refpos_cartesian( sub, anum2, : ) )
    end do
    end do

    !rotation along the xy axis
    !find the rotation angle
    rotax = (/0.d0,0.d0,1.d0/)
    rotangle = agg_angle
    !the lab frame
    call get_rmat( rmat, rotax, rotangle )
    to_cartesian = matmul( rmat, to_cartesian )
    !rotate the reference molecules
    do sub = 1, sublattice
    do anum2 = 1, anum
        refpos_cartesian( sub, anum2, : ) =    &
            matmul( rmat, refpos_cartesian( sub, anum2, : ) )
    end do
    end do

    !=================================================================!
    !build conversion matrices from frational to cartesian coordinates
    !=================================================================!
    
    !get the matrix to transform from fractional to cartesian
    call get_cartesian_mat( to_cartesian, cell_a, cell_b, cell_c,     &
                            cell_alpha, cell_beta, cell_gamma )
    !build the lattice along with a-axis
    n1 = 0
    do lb = 1, UC
    do sub = 1, sublattice
        n1 = n1 + 1
        dmat = (/0,lb-1,0/)
        mol1pos(n1,:,:) = refpos_cartesian( sub,:,:)
        mol1pos(n1,:,1) = mol1pos(n1,:,1) + dot_product( to_cartesian(1,:), dmat )
        mol1pos(n1,:,2) = mol1pos(n1,:,2) + dot_product( to_cartesian(2,:), dmat )
        mol1pos(n1,:,3) = mol1pos(n1,:,3) + dot_product( to_cartesian(3,:), dmat )
    end do
    end do
    ! print*, UC, sublattice
    ! print out the geometry
    ! do n1 = 1, nmax
    ! do da = 1, anum
        ! print*,mol1pos(n1,da,1),mol1pos(n1,da,2),mol1pos(n1,da,3)
    ! end do
    ! end do
    
end subroutine

!*********************************************************************!
!    Subroutine returns the matrix to transform to cartesian 
!    coordinates from fractional coordinates
!*********************************************************************!
subroutine get_cartesian_mat( to_cartesian, cell_a, cell_b, cell_c, &
                              cell_alpha, cell_beta, cell_gamma )
    implicit none
    
    real*8, intent(out) :: to_cartesian(3,3)
    real*8, intent(in) :: cell_a, cell_b, cell_c, cell_alpha,    &
                          cell_beta, cell_gamma
    
    to_cartesian = 0.d0
    to_cartesian(1,1) = cell_a
    to_cartesian(1,2) = cell_b * dcos( cell_gamma )
    to_cartesian(1,3) = cell_c * dcos( cell_beta  )
    to_cartesian(2,2) = cell_b * dsin( cell_gamma )
    to_cartesian(2,3) = cell_c * ( dcos( cell_alpha )   &
                                 - dcos( cell_beta )    &
                                 * dcos( cell_gamma ) ) &
                                 / dsin( cell_gamma )
    to_cartesian(3,3) = cell_c * dsqrt(                 &
                            1 - dcos(cell_alpha)**2     &
                              - dcos(cell_beta)**2      &
                              - dcos(cell_gamma)**2 +   &
                            2.d0*dcos(cell_alpha)       &
                                *dcos(cell_beta)        &
                                *dcos(cell_gamma) )     &
                                / dsin(cell_gamma)
end subroutine

subroutine get_rmat( rmat, rotax, rotangle )
    implicit none

    real*8, intent(out):: rmat(3,3)
    real*8, intent(in) :: rotax(3), rotangle
    real*8 :: kmat(3,3)
    
    !======================!
    !rodrigues formula
    !vrot = Rv
    !R = I + sin(q)K+(1-cos(q))K^2
    !v is the vector, I is identity matrix, R is rotation matrix, 
    !K is the cross product matrix of the rotation axis
    !            0    -kz    ky
    !    K =     kz    0    -kx
    !           -ky    kx    0
    !======================!
    kmat = 0.d0
    kmat(1,2) = -rotax(3)
    kmat(1,3) =  rotax(2)
    kmat(2,1) =  rotax(3)
    kmat(2,3) = -rotax(1)
    kmat(3,1) = -rotax(2)
    kmat(3,2) =  rotax(1)
    
    rmat = 0.d0
    rmat(1,1) = 1.d0
    rmat(2,2) = 1.d0
    rmat(3,3) = 1.d0
    rmat = rmat + dsin(rotangle)*kmat +             &
           (1-dcos(rotangle))*matmul(kmat,kmat)
end subroutine   

!*********************************************************************!
!    Subroutine returns the cross product of two 3 dimensional vectors
!*********************************************************************!
subroutine cross_product(vecina, vecinb, vecout)
    implicit none

    real*8, intent(in) :: vecina(3), vecinb(3)
    real*8, intent(out) :: vecout(3)
    
    vecout(1) = vecina(2)*vecinb(3) - vecina(3)*vecinb(2)
    vecout(2) = vecina(3)*vecinb(1) - vecina(1)*vecinb(3)
    vecout(3) = vecina(1)*vecinb(2) - vecina(2)*vecinb(1)
    
end subroutine