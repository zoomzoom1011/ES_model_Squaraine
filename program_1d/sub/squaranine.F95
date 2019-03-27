program squaraine
    use common_variables
    implicit none
            
    !read para file
    call read_para_file()
    ! call mkdir(trim(task_title))
    
    !index the basis states
    
    call index_mon_state()
    call index_sys_state()

    !allocate the hamiltonian eval
    call allocate_hev()
    
    !set the coupling and transition dipole moment
    call kount_out()
    call calc_coupling()
    
    !building the hamiltonian
    call hamiltonian()
        
    !output the hamiltonian 
    if (sys_kount < 500 )  call hamiltonian_out()
 
    !diagonalize
    call dev_diagonalize(sys_h, sys_kount, sys_eval)
   
    !output eigenvalues and eigenvectors
    if (sys_kount < 500 ) call evec_out()
        
    !calculate the absorption spectrum
    call charafxn()
    call abs_osc()
    call abs_spec()
    call abs_out()
    call normalization_ab()
        
    !calculate the pl spectrum
    if (calc_pl) call pl_osc()
    if (calc_pl) call pl_spec()
    if (calc_pl) call pl_out()
    if (calc_pl) call normalization_pl()

    !write the para to files
    call para_out()

end program