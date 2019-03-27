!*********************************************************************!
!       get the character of each state       !
!*********************************************************************!
subroutine charafxn()
    use common_variables
    implicit none
    
    integer:: fno
    character*256:: fname
    real*8:: fe(sys_kount)
    integer:: state, hx, n1, state_a
    
    fno = 101
    fname = trim(task_title)//'_chara.dat'  
    open( unit = fno, file = fname )
    
    !here i think the excite state mean one absorption
    fe = 0.d0
!    do state = 1, sys_kount
!        hx = 1, sys_kount
!        do n1 = 1, nmax
!            state_a = sys_state(hx,n1)%es_state
!            if (state_a == es_z1) then
!                fe(state) = fe(state) + sys_h(hx,state)**2
!            else if(state_a == es_z2) then
!                fe(state) = fe(state) + lambda_z2*sys_h(hx,state)**2
!            end if
!        end do
!        end do
        
        do state = 1, sys_kount
                fe(state) = fe(state) + sys_h(2,state)**2
                fe(state) = fe(state) + sys_h(3,state)**2
                fe(state) = fe(state) + sys_h(4,state)**2
                fe(state) = fe(state) + sys_h(7,state)**2
        write(fno, '(f14.7, 2f14.10)') (sys_eval(state)-sys_eval(1))* hw, fe(state)
    end do     
        
    close(fno)

end subroutine
