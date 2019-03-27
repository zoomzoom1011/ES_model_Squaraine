!*********************************************************************!
!        index system state                   !
!*********************************************************************!
subroutine  index_sys_state()
    use common_variables
    implicit none
    
#include<petsc/finclude/petscdef.h>
    integer:: run, n
    PetscInt:: sysnx
    PetscInt:: molecule_state(nmax)
    PetscInt, external:: get_mon_state, get_numvib
    real*8, external:: get_tot_charge

!    do i = 1, kount
!    do j = 1, nmax
!        mon_state = mod((i-1)/j, nmax) + 1
!        system_basis(i,j) = basis_set(mon_state)
!        print*, mon_state
!    end do
!    end do

!                state = 0
!        do m = 1, nmax
!            mon_state = 0
!    do estate = N, Z
!    do vib = 0, vibmax
!            state = state + 1
!                mon_state = mon_state + 1 
!                system_basis(state, m) = basis_set(mon_state)
!                print*, state, m,  mon_state
!        end do
!    end do
!    end do
    
    
    do run = 1, 2
        sys_kount = 0
        do sysnx = 1, mon_kount ** nmax
            do n = 1, nmax
                molecule_state(n) = get_mon_state(sysnx,n)
            end do
            
!            do n = 1, nmax
!                    print*, mon_state(molecule_state(n))%es_state
!                    print*, mon_state(molecule_state(n))%vib
!                        end do            
        
            if(get_numvib(molecule_state) > sys_vibmax) cycle
            ! if(get_tot_charge(molecule_state) .ne. 0) cycle
            sys_kount = sys_kount + 1
            
            if(run == 2) then
                do n = 1, nmax
                    sys_state(sys_kount,n) = mon_state(molecule_state(n))
!                    print*, sys_state(sys_kount,n)%es_state
!                    print*, sys_state(sys_kount,n)%vib
                end do
            end if
        end do
!                print*, sys_kount
        if(run == 1) then
            if(.not. allocated(sys_state)) allocate(sys_state(sys_kount, nmax))
        end if    
    end do
end subroutine

!*******************************************************************************!
!   return the molecule state for the given basis number & #molecule
!*******************************************************************************!
PetscInt function get_mon_state(sysnx,n)
    use common_variables
    implicit none
#include<petsc/finclude/petscdef.h>    
    PetscInt,intent(in):: sysnx
    integer,intent(in):: n
    
    get_mon_state = mod((sysnx-1)/(mon_kount**(n-1)),mon_kount)+1
    
    return
end function

!*******************************************************************************!
!   return the total molecule state vibrations for the given system basis 
!*******************************************************************************!
PetscInt function get_numvib(molecule_state)
    use common_variables
    implicit none
#include<petsc/finclude/petscdef.h>
    PetscInt,intent(in):: molecule_state(nmax)
    integer:: n
    
!    print*, molecule_state
    get_numvib = 0
    do n = 1, nmax
        get_numvib = get_numvib + mon_state(molecule_state(n))%vib1 + &
        mon_state(molecule_state(n))%vib2
    end do
    
!    print*, get_numvib
    return
end function

! !*******************************************************************************!
! !   return the molecule state for the given basis number & #molecule
! !*******************************************************************************!
! real*8 function get_tot_charge(molecule_state)
    ! use common_variables
    ! implicit none
! #include<petsc/finclude/petscdef.h>
    ! PetscInt, intent(in):: molecule_state(nmax)
    ! integer:: n, es_state, da
    
    ! get_tot_charge = 0.d0
    
    ! do n = 1, nmax
        ! es_state = mon_state(molecule_state(n))%es_state
        ! do da = 1, es_zmax
            ! get_tot_charge = get_tot_charge + tq(es_state, da)
        ! end do
    ! end do
    
! !    print*, get_tot_charge
    ! return
! end function
