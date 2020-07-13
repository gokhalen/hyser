include ../../src/HyserConf

# override implicit rules


${EQNOBJ}:${HYSER_COMMON_DEPS} ${EQNHDR} 
alleqn:${EQNOBJ}
	ar rcs libeqn.a ${EQNOBJ}

superclean:
	rm -rf *.o *.a


