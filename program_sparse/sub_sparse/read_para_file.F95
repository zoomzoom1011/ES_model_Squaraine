!*********************************************************************!
!                       read_para_file        !
!*********************************************************************!
subroutine read_para_file()
    use common_variables
    implicit none
    
    logical exists
    integer fno, ios, line, pos
    character*100 fname, buff, label
    
    fno = 201

    call getarg(1,fname)
    if (fname =="") then 
        if ( rank==0 ) print*, 'No control file given. Using default parameters'
        goto 1010
    end if 
    
!        make sure the control file exists
    inquire (file = trim (fname), exist = exists)       
    if (.not. exists) then 
        if ( rank==0 ) print*, 'Input file not found...aborting'
        stop
    end if
           
!        open the file to read in
    open (unit = fno, file= fname, status = 'old', action = 'read')
   
!--------------------------------------------------------------------------------------------------------! 
!        This process reads the control file
!        The ios changes if end of record or end of file is reached
!--------------------------------------------------------------------------------------------------------!

    ios = 0
    line = 0
    if ( rank==0 ) print*, 'Reading the input file...'
    
    do while(ios == 0)
!    print*, ios
        read(fno, '(a)', iostat = ios) buff
        if (ios == 0 ) then 
        line = line + 1 
        
        pos = scan(buff, ' ')
        label = buff(1:pos)
        buff = buff(pos + 1: )
        if (label (1:1)=='#' ) cycle
        select case (label)
            
        case('task_title')
            read(buff, *, iostat = ios) task_title
            if ( rank==0 ) print'(a)', 'setting task_title to:'//trim(task_title)
        case('es_zmax')
            read( buff, *, iostat=ios) es_zmax
            if ( rank==0 ) print'(a,i4)', 'Setting the number of ESM mode: ', es_zmax
        case( 'nmax' )
            read( buff, *, iostat=ios ) nmax
            if ( rank==0 ) print'(a,i4)', 'Setting nmax to: ', nmax
        case('vibmax')
            read( buff, *, iostat=ios ) vibmax
            if ( rank==0 ) print'(a,i4)', 'Setting vibmax to: ', vibmax
        case('sys_vibmax')
            read( buff, *, iostat=ios ) sys_vibmax
            if ( rank==0 ) print'(a,i4)', 'Setting sys_vibmax to: ', sys_vibmax
        case('hw')
            read( buff, *, iostat=ios ) hw
            if ( rank==0 ) print'(a,f10.2)', 'Setting vibration energy to: ', hw 
        ! case( 'sublattice' )
            ! read( buff, *, iostat=ios ) sublattice
            ! if ( rank==0 ) print'(a,i4)', 'Setting sublattice to: ', sublattice
        ! case('ct_state')
            ! read( buff, *, iostat=ios) ct_state
            ! if ( ct_state ) print'(a)', 'ct_states are turned on.'
            ! if ( .not.ct_state ) print'(a)', 'ct_states are turned off.'
!        case('s_z')
!            read( buff, *, iostat=ios ) s_z
!            print'(a,f8.2)', 'Setting s_z to : ', s_z
!        case('s_ani')
!            read( buff, *, iostat=ios ) s_ani
!            print'(a,f8.2)', 'Setting s_ani to : ', s_ani
!        case('s_cat')
!            read( buff, *, iostat=ios ) s_cat
!            print'(a,f8.2)', 'Setting s_cat to : ', s_cat    
        case('periodic')
            read(buff, *, iostat=ios) periodic
            if ( rank==0 )  then
                if (periodic) print'(a)', 'periodic condition is on.'
                if (.not. periodic) print'(a)', 'periodic condition is off.'
            end if
        case( 'lambda_n' )
            read( buff, *, iostat=ios ) lambda_n
            if ( rank==0 ) print'(a,f10.2)', 'Setting lambda_n to: ', lambda_n
        case( 'lambda_z1' )
            read( buff, *, iostat=ios ) lambda_z1
            if ( rank==0 ) print'(a,f10.2)', 'Setting lambda_z1 to: ', lambda_z1
        case( 'lambda_z2' )
            read( buff, *, iostat=ios ) lambda_z2
            if ( rank==0 ) print'(a,f10.2)', 'Setting lambda_z2 to: ', lambda_z2
        ! case( 'lambda_a' )
            ! read( buff, *, iostat=ios ) lambda_a
            ! if ( rank==0 ) print'(a,f10.2)', 'Setting lambda_a to: ', lambda_a
        ! case( 'lambda_c' )
            ! read( buff, *, iostat=ios ) lambda_c
            ! if ( rank==0 ) print'(a,f10.2)', 'Setting lambda_c to: ', lambda_c
        case('abs_lw')
            read( buff, *, iostat=ios) abs_lw
            if ( rank==0 ) print'(a,f8.2)', 'Setting the abs linewidth to (cm-1): ', abs_lw
        case('nz')
            read( buff, *, iostat=ios) nz
            if ( rank==0 ) print'(a,f8.2)', 'Setting zwitter energy to (cm-1): ', nz
        case('tz')
            read( buff, *, iostat=ios) tz
            if ( rank==0 ) print'(a,f8.2)', 'Setting intra charge transfer to (cm-1): ', tz
        ! case('nct')
            ! read( buff, *, iostat=ios) nct
            ! if ( rank==0 ) print'(a,f8.2)', 'Setting ion pair energy to (cm-1): ', nct
        ! case('tct')
            ! read( buff, *, iostat=ios) tct
            ! if ( rank==0 ) print'(a,f8.2)', 'Setting inter charge transfer to (cm-1): ', tct
        case('spec_step')
            read( buff, *, iostat=ios) spec_step
            if ( rank==0 ) print'(a,i8)', 'Setting spec_step to (cm-1): ', spec_step
        case('spec_start_ab')
            read( buff, *, iostat=ios) spec_start_ab
            if ( rank==0 ) print'(a,f8.2)', 'Setting spec_start_ab to (cm-1): ', spec_start_ab
        case('spec_end_ab')
            read( buff, *, iostat=ios) spec_end_ab
            if ( rank==0 ) print'(a,f8.2)', 'Setting spec_end_ab to (cm-1): ', spec_end_ab
        case('spec_start_pl')
            read( buff, *, iostat=ios) spec_start_pl
            if ( rank==0 ) print'(a,f8.2)', 'Setting spec_start_pl to (cm-1): ', spec_start_pl
        case('spec_end_pl')
            read( buff, *, iostat=ios) spec_end_pl
            if ( rank==0 ) print'(a,f8.2)', 'Setting spec_end_pl to (cm-1): ', spec_end_pl
        case('lorentzian')
            read(buff, *, iostat=ios) lorentzian
            if ( rank==0 )  then 
                if (lorentzian) print'(a)', 'Lineshape set to Lorentzian.'
                if (.not. lorentzian) print'(a)', 'Lineshape set to Gaussian.'
            end if
        case('calc_pl')
            read( buff, *, iostat=ios ) calc_pl
            if ( rank==0 ) then
                if ( calc_pl ) print'(a)', 'Will calculate all eigenstates.'
                if ( .not. calc_pl ) print'(a)', 'Will calculate only absorption.'
            end if
        case('x_offset')
            read( buff, *, iostat=ios) x_offset
            if ( rank==0 ) print'(a,f8.2)', 'Setting x_offset to (A): ', x_offset
        case('y_offset')
            read( buff, *, iostat=ios) y_offset
            if ( rank==0 ) print'(a,f8.2)', 'Setting y_offset to (A): ', y_offset
        case('z_offset')
            read( buff, *, iostat=ios) z_offset
            if ( rank==0 ) print'(a,f8.2)', 'Setting z_offset to (A): ', z_offset
        case('armlength')
            read( buff, *, iostat=ios) armlength
            if ( rank==0 ) print'(a,f8.2)', 'Setting armlength to (A): ', armlength
        case('agg_angle')
            read( buff, *, iostat=ios) agg_angle
            if ( rank==0 ) print'(a,f8.2)', 'Setting agg_angle to (degree): ', agg_angle
        case('dielectric')
            read( buff, *, iostat=ios) dielectric
            if ( rank==0 ) print'(a,f8.2)', 'Setting dielectric to (A): ', dielectric
        case( 'nev' )
            read( buff, *, iostat=ios ) nev
            if( rank==0 ) print'(a,i4)', 'Setting nev to: ', nev
        case('no_frenkel')
            read( buff, *, iostat=ios) no_frenkel
            if ( rank==0 )  then
                if ( no_frenkel ) print'(a)', 'no frenkel coupling will account.'
                if ( .not. no_frenkel ) print'(a)', 'frenkel coupling will account.'
            end if
        case('nearest_neighbor')
            read( buff, *, iostat=ios) nearest_neighbor
            if( rank==0 ) then
                if ( nearest_neighbor ) print'(a)', 'calc coupling from nearest neighbor.'
                if ( .not. nearest_neighbor ) print'(a)', 'calc coupling from long range.' 
            end if
        ! case('xyz_file1')
            ! read( buff, *, iostat=ios ) xyz_file(1)
            ! if( rank==0 ) print'(a)', 'Will read the xyz file: '//trim(xyz_file(1))
        ! case('xyz_file2')
            ! read( buff, *, iostat=ios ) xyz_file(2)
            ! if( rank==0 ) print'(a)', 'Will read the xyz file: '//trim(xyz_file(2))
        ! case('xyz_file3')
            ! read( buff, *, iostat=ios ) xyz_file(3)
            ! if( rank==0 ) print'(a)', 'Will read the xyz file: '//trim(xyz_file(3))
        ! case('xyz_file4')
            ! read( buff, *, iostat=ios ) xyz_file(4)
            ! if( rank==0 ) print'(a)', 'Will read the xyz file: '//trim(xyz_file(4))
        ! case('crystal_fname')
            ! read( buff, *, iostat=ios ) crystal_fname
            ! if( rank==0 ) print'(a)', 'Will read the crystal file: '//trim(crystal_fname)
        case default
            if ( rank==0 ) print'(a)', 'invalid label at line, ', line
            if ( rank==0 ) print'(a)', buff
            read*
        end select
        
        end if
        
    end do
    
!    print*, x_offset
!    print*, z_offset
    close(fno)

1010 continue
     if ( rank==0 ) print*, 'Calculating parameters in units of hw'

    !normalize parameters to hamiltonian units of vib quanta
    ! hw = hw*ev
    ! tz = tz*ev/hw
    ! nz = nz*ev/hw
    ! tct = tct*ev/hw
    ! nct = nct*ev/hw
    ! abs_lw = abs_lw*ev/hw
    ! spec_start_ab = spec_start_ab/hw
    ! spec_end_ab = spec_end_ab/hw
    ! spec_start_pl = spec_start_pl/hw
    ! spec_end_pl = spec_end_pl/hw
    
    ! if(ct_state) es_zmax = es_zmax + 2
    agg_angle = agg_angle/180.d0*pi
    
    !set lambda variables
    es_lambda_z1( es_n ) = lambda_n
    es_lambda_z1( es_z1 ) = lambda_z1
    es_lambda_z1( es_z2 ) = lambda_n
    es_lambda_z2( es_n ) = lambda_n
    es_lambda_z2( es_z1 ) = lambda_n
    es_lambda_z2( es_z2 ) = lambda_z2
    
    if ( vibmax == 0 ) then
        es_lambda_z1 = 0.d0
        es_lambda_z2 = 0.d0
    end if

    
end subroutine