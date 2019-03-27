!*********************************************************************!
!        calculate the absorption spectrum oscillator strength       !
!*********************************************************************!
subroutine abs_osc()
    use common_variables
    implicit none

    integer state, hx
    integer::fno = 202
        character*256::fname
    
    if(.not. allocated(ab_osc_x)) allocate (ab_osc_x(2:sys_kount))
    if(.not. allocated(ab_osc_y)) allocate (ab_osc_y(2:sys_kount))
    if(.not. allocated(ab_osc_z)) allocate (ab_osc_z(2:sys_kount))
    ab_osc_x = 0
    ab_osc_y = 0
    ab_osc_z = 0
        
        fname = trim(task_title)//'_abs_osc.dat'
        open(unit = fno, file = fname)
        
    do state = 2, sys_kount
        do hx = 1, sys_kount
            ab_osc_x(state) = ab_osc_x(state) + ux(hx) * sys_h(hx,state) * sys_h(hx,1)
            ab_osc_y(state) = ab_osc_y(state) + uy(hx) * sys_h(hx,state) * sys_h(hx,1)
            ab_osc_z(state) = ab_osc_z(state) + uz(hx) * sys_h(hx,state) * sys_h(hx,1)
        end do
        
        write(fno,'(i4,2f14.5)') state,(sys_eval(state)-sys_eval(1))*hw, &
            ab_osc_x(state)**2+ab_osc_z(state)**2
    end do
    
    
    ab_osc_x = ab_osc_x ** 2
    ab_osc_y = ab_osc_y ** 2
    ab_osc_z = ab_osc_z ** 2
end subroutine