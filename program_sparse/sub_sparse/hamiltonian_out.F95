!***********************************************************************!
!        output the hamiltonian to sparse matrix for calculation            !
!***********************************************************************!
subroutine hamiltonian_out()
    use common_variables
    implicit none
        
#include<petsc/finclude/petscdef.h>
    PetscInt state1, state2, dnz, onz, dnzmax, onzmax
    real*8, external :: get_matrix_element
    real*8 matrix_element, time_in_seconds
    integer run, diff
    integer, external:: get_diff
    
    if(rank == 0) print*, 'Building the sparse hamiltonian matrix' 
    time_in_seconds = MPI_Wtime()
    
    !create the sparse matrix by communicator
    call MatCreate(PETSC_COMM_WORLD, h_sparse, ierr)
    !Sets the local and global sizes, and checks to determine compatibility
    call MatSetSizes(h_sparse, PETSC_DECIDE, PETSC_DECIDE, sys_kount, sys_kount, ierr)
    !Creates a matrix where the type is determined from the options database
    call MatSetFromOptions(h_sparse, ierr)
    !Sets up the internal matrix data structures for the later use.
    call MatSetUp(h_sparse, ierr)
    !Returns the range of matrix rows owned by this processor
    call MatGetOwnershipRange(h_sparse,Istart,Iend,ierr)
    !make sure it's symmetric in terms of both structure and value
    call MatSetOption(h_sparse, MAT_SYMMETRIC, PETSC_TRUE, ierr )
    
     !initialize
     dnzmax = 0
     onzmax = 0

    !the first run, you preallocate the matrix, the second you actually
    !build it
    do run = 1, 2
        do state1 = Istart + 1, Iend
            !only create on the local process
            if ( rank == 0 ) then
                !update user on progress
                print'(1a1,i4,a2,i4,a1,$)', char(13), run,'/2', 100*state1/(Iend-Istart), '%'
            end if
            
            !counter for nonzero's, initialize
            dnz = 0
            onz = 0
            do state2 = 1, sys_kount
                diff = get_diff( state1, state2 )
                if ( diff > 1 ) cycle !try to speed this up by only calculating matrix
                                      !elements for states you no won't be zero
                matrix_element = get_matrix_element( state1, state2 )
                if ( matrix_element .ne. 0.d0 ) then
                    !count nonzero elements for matrix allocation and filling
                    if ( state2 > Istart .and. state2<= Iend ) then
                        dnz = dnz + 1
                    else
                        onz = onz + 1
                    end if
                    !set the matrix value...I was having problems with MatSetValues, so
                    !I just will do one at a time
                    if ( run == 2 ) then
                        call MatSetValue(h_sparse, state1-1, state2-1, matrix_element, &
                                          INSERT_VALUES, ierr )
                    end if
                end if
            end do

            if ( run == 1 ) then
            !set the total diagonal and off diagonal nonzeros for each row
            !and allocate arrays for setting matrix values
                if ( dnz > dnzmax ) dnzmax = dnz
                if ( onz > onzmax ) onzmax = onz
            end if

        end do
        
        !preallocate the matrices
        !set the total diagonal and off diagonal nonzeros for each row
        !and allocate arrays for setting matrix values
        if ( run == 1 ) then
            if ( rankmax == 1 ) then
                call MatSeqAIJSetPreallocation(h_sparse, dnzmax, PETSC_NULL_INTEGER, ierr )
            else
                call MatMPIAIJSetPreallocation(h_sparse, dnzmax, PETSC_NULL_INTEGER, &
                                                onzmax, PETSC_NULL_INTEGER, ierr )
            end if
        end if
    end do
    
    call MatAssemblyBegin( h_sparse, MAT_FINAL_ASSEMBLY, ierr )
    call MatAssemblyEnd( h_sparse, MAT_FINAL_ASSEMBLY, ierr )

    time_in_seconds = MPI_Wtime() - time_in_seconds
    if (rank==0) print'(a,f10.4)', 'Done... Time: ', time_in_seconds

end subroutine

!**********************************************************************!
integer function get_diff( state1, state2 )
    use common_variables
    implicit none

#include<petsc/finclude/petscdef.h>
    PetscInt, intent(in) :: state1, state2
    integer state_a, vib_a1, vib_a2, &
            state_b, vib_b1, vib_b2, &
            n

    get_diff = 0
    do n = 1, nmax
        state_a = sys_state(state1, n)%es_state
        vib_a1 = sys_state(state1, n)%vib1
        vib_a2 = sys_state(state1, n)%vib2
        state_b = sys_state(state2, n)%es_state
        vib_b1 = sys_state(state2, n)%vib1
        vib_b2 = sys_state(state2, n)%vib2

        !check if the state of the current molecule is different in a
        !and b
        if ( vib_a1 .ne. vib_b1 .or. & 
             vib_a2 .ne. vib_b2 .or. & 
             state_a .ne. state_b ) then
             get_diff = get_diff + 1
             !at max, only one molecule can have a different
             !configuration
             if ( get_diff > 1 ) return
        end if
    end do
end function

                
!*********************************************************************!
!         build hamiltonian fucnction for sparse matrix technique                   !
!*********************************************************************!
real*8 function get_matrix_element(state1, state2)
    use common_variables
    implicit none

#include<petsc/finclude/petscdef.h>
    PetscInt, intent(in) :: state1, state2
    integer:: n1
    integer:: state_a, vib_a1, vib_a2
    integer:: state_b, vib_b1, vib_b2
    integer:: diff, diffn
    integer:: diffterm(3)
    real*8:: volap, volapfact
    integer:: state_1a, state_2a, state_1b, state_2b
    integer:: vib_1a, vib_2a, vib_1b, vib_2b
    ! real*8:: get_tct
    
    get_matrix_element = 0.d0
    
    !on diagonal
    if(state1 == state2) then
        do n1 = 1, nmax
            state_a = sys_state(state1,n1)%es_state
            vib_a1 = sys_state(state1,n1)%vib1
            vib_a2 = sys_state(state1,n1)%vib2
            
            !the electronic energy
            if (state_a == es_n) then
                get_matrix_element = get_matrix_element + 0.d0
            else if (state_a == es_z1) then
                get_matrix_element = get_matrix_element + nz
            else if (state_a == es_z2) then
                get_matrix_element = get_matrix_element + nz
            end if

            !the vibration energy
            get_matrix_element = get_matrix_element + (vib_a1+vib_a2)*hw
        end do
    
        !coupling for the dimer
        if(nmax > 1) then
            get_matrix_element =  get_matrix_element + coulomb_coupling(state1)
        end if
    
    !off diagonal
    else
        diff = 0
        do n1 = 1, nmax
            state_a = sys_state(state1,n1)%es_state
            vib_a1 = sys_state(state1,n1)%vib1
            vib_a2 = sys_state(state1,n1)%vib2
            state_b = sys_state(state2,n1)%es_state
            vib_b1 = sys_state(state2,n1)%vib1
            vib_b2 = sys_state(state2,n1)%vib2
                    
            if ( vib_a1 .ne. vib_b1 .or. & 
                 vib_a2 .ne. vib_b2 .or. & 
                 state_a .ne. state_b ) then
            
                diffn = n1
                diff = diff + 1
                !at max, only one molecule can have a different
                !configuration
                if ( diff > 1 ) return
            end if
        end do
        
        !at max, only one molecule can have a different configuration
        if ( diff == 1 ) then
            
            state_a = sys_state(state1,diffn)%es_state
            vib_a1 = sys_state(state1,diffn)%vib1
            vib_a2 = sys_state(state1,diffn)%vib2
            state_b = sys_state(state2,diffn)%es_state
            vib_b1 = sys_state(state2,diffn)%vib1
            vib_b2 = sys_state(state2,diffn)%vib2
            
            !intramolecule coupling
            ! if ( state_a == es_n .and. state_b == es_z) then
                ! get_matrix_element = get_matrix_element  + tz* volap( es_lambda_z( state_a ), vib_a,&
                                       ! es_lambda_z( state_b ), vib_b)
            ! else if ( state_a == es_z .and. state_b == es_n ) then
                ! get_matrix_element = get_matrix_element + tz* volap( es_lambda_z( state_a ), vib_a,&
                                       ! es_lambda_z( state_b ), vib_b)
            ! end if

            if ( state_a == es_n .and. state_b == es_z1 ) then
                get_matrix_element = tz
            else if ( state_a == es_n .and. state_b == es_z2 ) then
                get_matrix_element = tz
            else if ( state_a == es_z1 .and. state_b == es_n ) then
                get_matrix_element = tz
            else if ( state_a == es_z2 .and. state_b == es_n ) then
                get_matrix_element = tz
            else
                get_matrix_element = 0.d0
                ! return
            end if
            !put in the vibrations last since these are bigger calcs
            volapfact = volap( es_lambda_z1( state_a ), vib_a1,&
                               es_lambda_z1( state_b ), vib_b1)&
                       *volap( es_lambda_z2( state_a ), vib_a2,&
                               es_lambda_z2( state_b ), vib_b2)
            ! if(volapfact .eq. 0.d0) return
            get_matrix_element = get_matrix_element*volapfact
            return
        end if 
        
    !intermolecule coupling

    ! if(nmax > 1 and ct_state) then

        ! diff = 0
        ! do n1 = 1, nmax
            ! state_a = sys_state(state1,n1)%es_state
            ! vib_a = sys_state(state1,n1)%vib
            ! state_b = sys_state(state2,n1)%es_state
            ! vib_b = sys_state(state2,n1)%vib
                    
            ! if ( vib_a .ne. vib_b .or. & 
                 ! state_a .ne. state_b ) then
                
                ! diff = diff + 1
                ! diffterm(diff) = n1
            ! end if
            ! if ( diff > 2 ) cycle

        ! end do

        ! !at max, only two molecules can have different configurations

        ! if ( diff == 2 ) then !make sure it is charge transfer state

            ! state_1a = sys_state(state1,diffterm(1))%es_state
            ! vib_1a = sys_state(state1,diffterm(1))%vib
            ! state_1b = sys_state(state1,diffterm(2))%es_state
            ! vib_1b = sys_state(state1,diffterm(2))%vib
            ! state_2a = sys_state(state2,diffterm(1))%es_state
            ! vib_2a = sys_state(state2,diffterm(1))%vib
            ! state_2b = sys_state(state2,diffterm(2))%es_state
            ! vib_2b = sys_state(state2,diffterm(2))%vib
            ! diff = 0
            ! if (state_1a == es_a) diff = diff + 1
            ! if (state_1b == es_a) diff = diff + 1
            ! if (state_2a == es_a) diff = diff + 1
            ! if (state_2b == es_a) diff = diff + 1
            ! if(diff == 1) then
                ! !vibronic coupling

                ! get_matrix_element = get_matrix_element + get_tct(diffterm(1),diffterm(2)) * &
                                            ! volap( es_lambda_z( state_1a ), vib_1a,&
                                            ! es_lambda_z( state_2a ), vib_2a) * &
                                            ! volap( es_lambda_z( state_1b ), vib_1b,&
                                            ! es_lambda_z( state_2b ), vib_2b)
                                        
        ! !        get_matrix_element = get_matrix_element + tct* &
        ! !                                    volap( es_lambda_z( state_1a ), vib_1a,&
        ! !                                    es_lambda_z( state_2a ), vib_2a) * &
        ! !                                    volap( es_lambda_z( state_1b ), vib_1b,&
        ! !                                    es_lambda_z( state_2b ), vib_2b)
            ! end if
        ! end if
    ! end if
        
    
    end if
end function