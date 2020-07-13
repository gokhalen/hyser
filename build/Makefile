include ./HyserConf


HYSERLIBS=${ELEMLIB} ${INTELIB} ${SHAPLIB} ${UTILLIB} ${DATLIB} ${MESHLIB} ${CONSLIB} ${EQNLIB}

${HYSERLIBS}:${LIBSRC} ${LIBHDR} ${HYSER_COMMON_DEPS}
	@cd ${LIBDIR}; ${MAKE} alllib

${FORWARDOBJ}:${LIBSRC} ${LIBHDR} ${HYSER_COMMON_DEPS}
${INVERSEOBJ}:${LIBSRC} ${LIBHDR} ${HYSER_COMMON_DEPS}

forward: ${FORWARDOBJ} ${HYSERLIBS} ${HYSER_COMMON_DEPS}
	@echo ${FORWARDOBJ}
	${CLINKER} ${SLEPC_OPT}  -o forward ${FORWARDOBJ} ${CPPFLAGS}$ ${PETSC_SNES_LIB} ${SLEPC_LIB} ${HYSERLIBS}

#inverse: ${INVERSEOBJ} ${HYSERLIBS} ${HYSER_COMMON_DEPS}
#	@echo ${INVERSEOBJ}
#	${CLINKER} ${SLEPC_OPT}  -o inverse ${INVERSEOBJ} ${CPPFLAGS}$ ${PETSC_SNES_LIB} ${SLEPC_LIB} ${TAO_LIB} ${HYSERLIBS}


superclean:
#	rm - rf forward inverse
	@echo ${HYSER_DIR}
	@cd ${LIBDIR};  ${MAKE} superclean; cd ${SRCDIR}
	rm -rf *.o *.exe *~* forward 
# DO NOT DELETE

printvars:
	@echo "CLINKER= ${CLINKER}"

tar:
	@cd ${SRCDIR};cd ..; rm hyserIII.tar.gz; tar -czvf hyserIII.tar.gz docs headers/ lib/ src/ test/ utils/ ;cd ${SRCDIR}

