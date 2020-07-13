#Include top level makefile. This gives flags and variables
include ../../src/HyserConf
${INTEOBJ}:${HYSER_COMMON_DEPS} ${INTEHDR}


allint:${INTEOBJ}
	ar rcs libinteg.a ${INTEOBJ}

superclean:
	rm -rf *.o *.a
