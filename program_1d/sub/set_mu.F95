!*********************************************************************!
!        set transition dipole moment to state !
!*********************************************************************!
subroutine set_mu()
    use common_variables
    implicit none
    
    integer i, j
    
    if(.not. allocated(ux)) allocate(ux(sys_kount))
    if(.not. allocated(uy)) allocate(uy(sys_kount))
    if(.not. allocated(uz)) allocate(uz(sys_kount))
    ux = 0.d0
    uy = 0.d0
    uz = 0.d0
        
    do i = 1, sys_kount
    do j = 1, nmax
        if (sys_state(i,j)%es_state == 2) then
            ux(i) = ux(i) + 1.d0
        else if (sys_state(i,j)%es_state == 3) then
            ux(i) = ux(i) - 1.d0
        end if
    end do
    end do
    
end subroutine

