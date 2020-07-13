#Include top level makefile. This gives flags and variables
include ../../src/HyserConf


${UTILOBJ}:${HYSER_COMMON_DEPS} ${UTILHDR}
allutil:${UTILOBJ}
	ar rcs libutil.a ${UTILOBJ}

superclean:
	rm -rf *.a *.o


