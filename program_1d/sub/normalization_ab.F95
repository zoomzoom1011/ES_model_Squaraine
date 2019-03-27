!*********************************************************************!
!    Normalization for the abs spectrum
!*********************************************************************!
subroutine normalization_ab()
        use common_variables
        implicit none
        
        real*8 maxvalue
        character*256 f_name
    real*8 energy
    integer spec_point, f_no
                
        maxvalue = 0
        do spec_point=1, spec_step
                if (ab_x(spec_point) > maxvalue) maxvalue=ab_x(spec_point)
        end do
        
        ab_x = ab_x / maxvalue
    
    f_no = 205
    f_name = trim(task_title)//'_ab_normal.dat'
    
    open(unit=f_no, file=f_name)
        write( f_no, * ) 'Energy A(\g(w))'
                write( f_no, * ) 'Energy System'
                write( f_no, * ) 'cm\+(-1) a.u.'
        do spec_point = 1, spec_step
            energy = spec_start_ab + (spec_end_ab - spec_start_ab)/(spec_step*1.d0)*spec_point
            write(f_no,'(14f14.7)') energy*hw, ab_x(spec_point)
        end do
    
    close(f_no)
                          
end subroutine