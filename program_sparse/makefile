myOBJDIR=obj
mySRCDIR=sub_sparse
mySRCFILES=squaraine.F95 abs_out.F95 abs_osc.F95 abs_spec.F95 \
		hamiltonian_out.F95 bring_inside_aggregate_range.F95 \
		index_mon_state.F95 index_sys_state.F95 \
		kount_out.F95 \
		normalization_ab.F95 \
		para_out.F95 \
		read_para_file.F95 volap.F95 \
		normalization_pl.F95 normalization_ab.F95 \
		build_molecule_geometry.F95 calc_coupling.F95 \
		pl_osc.F95 pl_out.F95 pl_spec.F95
myOBJFILES=$(mySRCFILES:.F95=.o)
myEXE=squaraine.exe
myDEP=common_variables.mod
vpath %.o ${myOBJDIR}
vpath %.mod ${myOBJDIR}
vpath %.F95 ${mySRCDIR}
FLINKER=gfortran
SLEPC_DIR=${HOME}/slepc
PETSC_DIR=${HOME}/petsc
include ${SLEPC_DIR}/lib/slepc/conf/slepc_variables


${myEXE} : ${myOBJFILES} ${myDEP:.mod=.o}
	-${FLINKER} $(addprefix ${myOBJDIR}/,$(^F)) -o $@ -I${myOBJDIR} ${SLEPC_EPS_LIB}
	# -${FLINKER} $(addprefix ${myOBJDIR}/,$(^F)) -o $@ -I${myOBJDIR}
%.o : %.F95 ${myDEP:.mod=.o}
	-${FLINKER} -c ${mySRCDIR}/$(<F) -o ${myOBJDIR}/$(@F) -I${myOBJDIR} ${PETSC_FC_INCLUDES} ${SLEPC_FC_INCLUDES}
	# -${FLINKER} -c ${mySRCDIR}/$(<F) -o ${myOBJDIR}/$(@F) -I${myOBJDIR}
${myDEP:.mod=.o} : ${myDEP:.mod=.F95}
	-${FLINKER} -c ${mySRCDIR}/${myDEP:.mod=.F95} -o ${myOBJDIR}/${myDEP:.mod=.o} -J${myOBJDIR} ${PETSC_FC_INCLUDES} ${SLEPC_FC_INCLUDES}
	# -${FLINKER} -c ${mySRCDIR}/${myDEP:.mod=.F95} -o ${myOBJDIR}/${myDEP:.mod=.o} -J${myOBJDIR}
clean:
	@echo cleaning up
	@-rm ${myOBJDIR}/*.o ${myOBJDIR}/*.mod 2>/dev/null || true
	@-rm ${myEXE}
	@clear
