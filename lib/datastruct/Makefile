#Include top level makefile. This gives flags and variables
include ../../src/HyserConf

${DATOBJ}:${HYSER_COMMON_DEPS} ${DATHDR}

alldat:${DATOBJ} 
	ar rcs libdat.a ${DATOBJ}

superclean:
	rm -rf *.o *.a
