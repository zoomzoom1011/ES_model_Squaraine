!*********************************************************************!
!         build hamiltonian                    !
!*********************************************************************!
subroutine hamiltonian()
    use common_variables
    implicit none
    
    integer:: state1,state2
    integer:: n1
    integer:: state_a, vib_a1, vib_a2
    integer:: state_b, vib_b1, vib_b2
    integer:: diff, diffn
    real*8:: volap, volapfact
    
    do state1 = 1, sys_kount
    
        !on diagonal
        do n1 = 1, nmax
                        state_a = sys_state(state1,n1)%es_state
            vib_a1 = sys_state(state1,n1)%vib1
                        vib_a2 = sys_state(state1,n1)%vib2
            
            !the electronic energy
            if (state_a == es_z1) then
                sys_h(state1,state1)=sys_h(state1,state1) + etaz
            else if (state_a == es_z2) then
                sys_h(state1,state1)=sys_h(state1,state1) + etaz
            end if
            
            !the vibration energy
            sys_h(state1, state1) = sys_h(state1, state1) + (vib_a1 + vib_a2)*1.d0
        end do
        
        !coupling for the dimer
        if(nmax > 1) then
            sys_h(state1, state1) =  sys_h(state1, state1) + coulomb_coupling(state1)
        end if
    end do
    
    !off diagonal
    do state1 = 1 ,sys_kount
    do state2 = 1 ,sys_kount
        if(state1 == state2) cycle
        !intramolecular CT
        diff = 0
        do n1 = 1, nmax
            state_a = sys_state(state1,n1)%es_state
            vib_a1 = sys_state(state1,n1)%vib1
                        vib_a2 = sys_state(state1,n1)%vib2
            state_b = sys_state(state2,n1)%es_state
            vib_b1 = sys_state(state2,n1)%vib1
                        vib_b2 = sys_state(state2,n1)%vib2
                    
            if ( vib_a1 .ne. vib_b1 .or. & 
                 state_a .ne. state_b .or. & 
                                 vib_a2 .ne. vib_b2 ) then
            
                diffn = n1
                diff = diff + 1
            end if
            if ( diff > 1 ) cycle
        end do
        
        !at max, only one molecule can have a different configuration
        if ( diff .ne. 1 ) cycle
        
        state_a = sys_state(state1,diffn)%es_state
        vib_a1 = sys_state(state1,diffn)%vib1
        vib_a2 = sys_state(state1,diffn)%vib2
        state_b = sys_state(state2,diffn)%es_state
        vib_b1 = sys_state(state2,diffn)%vib1
        vib_b2 = sys_state(state2,diffn)%vib2
        
        
        if ( state_a == es_n .and. state_b == es_z1 ) then
            sys_h(state1, state2) = tz
        else if ( state_a == es_n .and. state_b == es_z2 ) then
            sys_h(state1, state2) = tz
        else if ( state_a == es_z1 .and. state_b == es_n ) then
            sys_h(state1, state2) = tz
        else if ( state_a == es_z2 .and. state_b == es_n ) then
            sys_h(state1, state2) = tz
        else
            sys_h(state1, state2) = 0.d0
        end if
        
        !put vibration into
        volapfact = volap( es_lambda_z1( state_a ), vib_a1,&
                           es_lambda_z1( state_b ), vib_b1)&
                   *volap( es_lambda_z2( state_a ), vib_a2,&
                           es_lambda_z2( state_b ), vib_b2)
        sys_h(state1, state2) = sys_h(state1, state2) * volapfact
        ! print*, state_a, vib_a1, vib_a2, state_b, vib_b1, vib_b2
    end do
    end do
            

end subroutine