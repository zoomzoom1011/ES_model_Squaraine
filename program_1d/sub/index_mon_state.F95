!*********************************************************************!
!        index the basis state      !
!*********************************************************************!
subroutine  index_mon_state()
    use common_variables
    implicit none
    
    integer es_state
    integer vib1, vib2
    integer run
    
!    kount_state = 0
!    do estate = 1, es_zmax
!    do vib = 0, vibmax
!        kount_state = kount_state + 1
!    end do 
!    end do
    
!    state = 0
!    if (.not. allocated (basis_set)) allocate(basis_set(kount_state)) 
    
!    do estate = 1, es_zmax
!    do vib = 0, vibmax
!        state = state + 1
!        basis_set(state)%electronic_state = estate
!        basis_set(state)%vib_state = vib
!        print*, vib, estate
!    end do 
!    end do

    do run = 1, 2
    mon_kount = 0
        do es_state = 1, es_zmax
        do vib1 = 0, vibmax
            do vib2 = 0, vibmax
            mon_kount = mon_kount + 1
            if(run == 2) then
                mon_state(mon_kount)%es_state = es_state
                mon_state(mon_kount)%vib1 = vib1
                mon_state(mon_kount)%vib2 = vib2
!                    print*, es_state, vib
            end if
        end do 
        end do
        end do
        
    if (run == 1) then
        if (.not. allocated (mon_state)) allocate(mon_state(mon_kount)) 
    end if
    end do
    
end subroutine