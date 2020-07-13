#Include top level Conf. This gives flags and variables
include ../../src/HyserConf

DEPS=../datastruct/*.cpp 
${ELEMOBJ}:${HYSER_COMMON_DEPS} ${DEPS}  ${ELEMHDR}

allelm:${ELEMOBJ} 
	ar rcs libelem.a ${ELEMOBJ}

superclean:
	${RM} -rf *.a *.o

