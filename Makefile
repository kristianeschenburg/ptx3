include ${FSLCONFDIR}/default.mk

PROJNAME = fdt

# Commt this line out to compile without using debugger
OPTFLAGS = -ggdb

#ARCHFLAGS = -arch i386
#ARCHLDFLAGS = -arch i386

ifeq ($(FSLMACHTYPE),apple-darwin8-gcc4.0)
        ARCHFLAGS =  -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -I/usr/X11R6/include/
        ARCHLDFLAGS = -Wl,-search_paths_first -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -L/Developer/SDKs/MacOSX10.4u.sdk/usr/X11R6/lib/ 
endif 


USRINCFLAGS = -I ${INC_EXPAT} -I ${INC_NEWMAT} -I ${INC_NEWRAN} -I ${INC_CPROB} -I ${INC_PROB} -I ${INC_BOOST} -I ${INC_ZLIB} -I /mnt/parcellator/parcellation/Software/fslbuild/fsl/src/
USRLDFLAGS = -L ${LIB_EXPAT} -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

DLIBS = -lnewmeshclass -lwarpfns -lbasisfield -lfslsurface  -lfslvtkio -lmeshclass -lnewimage -lutils -lmiscmaths -lnewmat -lnewran -lfslio -lgiftiio -lexpat -lfirst_lib -lniftiio -lznz -lutils -lprob -lm -lz

LINK = -L/mnt/parcellator/parcellation/Software/fslbuild/fsl/src/ptx3/

CCOPS=ccops
PTX=probtrackx3
FTB=find_the_biggest
PJ=proj_thresh
FMO=fdt_matrix_ops
FMS=fdt_matrix_split
TEST=testfile

CCOPSOBJS=ccops.o ccopsOptions.o 
PTXOBJS=probtrackx.o probtrackxOptions.o streamlines.o ptx_simple.o ptx_seedmask.o ptx_nmasks.o csv.o csv_mesh.o
FTBOBJS=find_the_biggest.o csv_mesh.o
PJOBJS=proj_thresh.o csv_mesh.o
FMOOBJS=fdt_matrix_ops.o
FMSOBJS=fdt_matrix_split.o
TESTOBJS=testfile.o streamlines.o csv.o csv_mesh.o probtrackxOptions.o

SURFDATA=surf_proj
SURFDATAOBJS=surf_proj.o csv.o csv_mesh.o
MATMERGE=fdt_matrix_merge
MATMERGEOBJS=fdt_matrix_merge.o streamlines.o csv.o csv_mesh.o probtrackxOptions.o
MAT42=fdt_matrix_4_to_2
MAT42OBJS=fdt_matrix_4_to_2.o 

S2S=surf2surf
S2SOBJS=surf2surf.o csv.o csv_mesh.o
S2V=surf2volume
S2VOBJS=surf2volume.o csv.o csv_mesh.o
L2S=label2surf
L2SOBJS=label2surf.o csv_mesh.o
SM=surfmaths
SMOBJS=surfmaths.o 


XFILES = probtrackx3


all: ${XFILES} 

${CCOPS}:    	${CCOPSOBJS}	
		   ${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${CCOPSOBJS} ${DLIBS}

${FTB}:    	${FTBOBJS}
		   ${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${FTBOBJS} ${DLIBS} 

${PJ}:    	${PJOBJS}
		   ${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PJOBJS} ${DLIBS} 


${FMO}:    	${FMOOBJS}
		   ${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${FMOOBJS} ${DLIBS}

${FMS}:    	${FMSOBJS}
		   ${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${FMSOBJS} ${DLIBS}

${TEST}:    	${TESTOBJS}
		   ${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${TESTOBJS} ${DLIBS}


${SURFDATA}:    ${SURFDATAOBJS}
		${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SURFDATAOBJS} ${DLIBS}

${MATMERGE}:    ${MATMERGEOBJS}
		${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MATMERGEOBJS} ${DLIBS}

${MAT42}:    	${MAT42OBJS}
		${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MAT42OBJS} ${DLIBS}


${S2S}:         ${S2SOBJS}
		${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${S2SOBJS} ${DLIBS}

${S2V}:         ${S2VOBJS}
		${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${S2VOBJS} ${DLIBS}

${L2S}:         ${L2SOBJS}
		${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${L2SOBJS} ${DLIBS}

${SM}:         ${SMOBJS}
		${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SMOBJS} ${DLIBS}
		
${PTX}:		   ${PTXOBJS}
		   ${CXX} ${LINK} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PTXOBJS} ${DLIBS}

