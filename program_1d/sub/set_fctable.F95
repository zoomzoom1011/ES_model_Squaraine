subroutine set_fctable()
    use common_variables
    implicit none
    
    real*8 fc, s_diff
    integer vibg,vib,vibc,viba
    
    if ( .not. allocated(fc_gz)) allocate (fc_gz(0:vibmax, 0:vibmax))
    if ( .not. allocated(fc_gc)) allocate (fc_gc(0:vibmax, 0:vibmax))
    if ( .not. allocated(fc_ga)) allocate (fc_ga(0:vibmax, 0:vibmax))
    if ( .not. allocated(fc_cz)) allocate (fc_cz(0:vibmax, 0:vibmax))
    if ( .not. allocated(fc_az)) allocate (fc_az(0:vibmax, 0:vibmax))
    fc_gz = 0.d0
    fc_gc = 0.d0
    fc_ga = 0.d0
    fc_cz = 0.d0
    fc_az = 0.d0
    
    !ground to excited state
    do vibg = 0 , vibmax
    do vib = 0,  vibmax
        call fcfac(vibg, vib, s_z, fc)
        fc_gz(vibg,vib) = fc
    end do
    end do
    
    !ground to cation state
    do vibg = 0 , vibmax
    do vibc = 0,  vibmax
        call fcfac(vibg, vibc, s_cat, fc)
        fc_gc(vibg,vibc) = fc
    end do
    end do        

    !ground to anion state
    do vibg = 0, vibmax
    do viba = 0, vibmax
        call fcfac(vibg, viba, s_ani, fc)
        fc_ga(vibg,viba) = fc
    end do
    end do
    
    !cation to zwitter    
    do vibc = 0,vibmax
    do vib  = 0,vibmax
        s_diff = dsqrt(s_z) - dsqrt(s_cat)
        s_diff = s_diff**2
        if( s_z <= s_cat ) then
            call fcfac(vib,vibc,dabs( s_diff ),fc)
        else
            call fcfac(vibc,vib,dabs( s_diff ),fc)
        endif
        
        fc_cz(vibc, vib) = fc
    enddo
    enddo
    
    !anion to zwitter
    do viba = 0,vibmax
    do vib = 0,vibmax
        s_diff = dsqrt(s_z) - dsqrt(s_ani)
        s_diff = s_diff**2
        if( s_z <= s_ani ) then
            call fcfac(vib,viba,dabs( s_diff ),fc)
        else
            call fcfac(viba,vib,dabs( s_diff ),fc)
        endif
        fc_az(viba, vib) = fc
    enddo
    enddo

end subroutine