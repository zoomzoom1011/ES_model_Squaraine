!*********************************************************************!
!                pl spectrum to each state         !
!*********************************************************************!
subroutine pl_osc()
    use common_variables
    implicit none
    
    integer hx
    integer state, state1
    integer::fno = 202
        character*256::fname
        
    if(.not. allocated(pl_osc_x)) allocate(pl_osc_x(sys_kount, sys_kount))
    if(.not. allocated(pl_osc_y)) allocate(pl_osc_y(sys_kount, sys_kount))
    if(.not. allocated(pl_osc_z)) allocate(pl_osc_z(sys_kount, sys_kount))
    pl_osc_x = 0.d0
    pl_osc_y = 0.d0
    pl_osc_z = 0.d0
    
        fname = trim(task_title)//'_pl_osc.dat'
        open(unit = fno, file = fname)
        
    do state = 2, sys_kount
    do state1 = 1, state-1
            do hx = 1, sys_kount
                pl_osc_x(state1,state) = pl_osc_x(state1,state) + ux(hx) * sys_h(hx,state) * sys_h(hx,state1)
                pl_osc_y(state1,state) = pl_osc_y(state1,state) + uy(hx) * sys_h(hx,state) * sys_h(hx,state1)
                pl_osc_z(state1,state) = pl_osc_z(state1,state) + uz(hx) * sys_h(hx,state) * sys_h(hx,state1)
            end do
    
        write(fno,'(2i4,2f14.6)') state,state1,(sys_eval(state)-sys_eval(state1))*hw, &
                pl_osc_x(state1,state)**2+pl_osc_z(state1,state)**2
        end do
    end do

    pl_osc_x = pl_osc_x ** 2
    pl_osc_y = pl_osc_y ** 2
    pl_osc_z = pl_osc_z ** 2
end subroutine