include ../src/HyserConf

alllib: 
	@cd ${UTILDIR}; ${MAKE} allutil
	@cd ${DATDIR};  ${MAKE} alldat;
	@cd ${ELEMDIR}; ${MAKE} allelm;
	@cd ${INTEDIR}; ${MAKE} allint;
	@cd ${SHAPDIR}; ${MAKE} allshp;
	@cd ${MESHDIR}; ${MAKE} allmsh;
	@cd ${EQNDIR};	${MAKE} alleqn;
	@cd ${CONSDIR};	${MAKE} allcons;
	@cd ${HYSER_DIR}/lib;              

superclean:
	@cd ${UTILDIR}; ${MAKE} superclean;
	@cd ${DATDIR};  ${MAKE} superclean;   
	@cd ${ELEMDIR}; ${MAKE} superclean;    
	@cd ${INTEDIR}; ${MAKE} superclean;    
	@cd ${SHAPDIR}; ${MAKE} superclean;    
	@cd ${MESHDIR}; ${MAKE} superclean;    
	@cd ${EQNDIR};	${MAKE} superclean;    
	@cd ${CONSDIR};	${MAKE} superclean;    
	@cd ${HYSER_DIR}/lib;              