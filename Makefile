all:
	${MAKE} -f psar.mk
	cd fsa-1.15.6 && ./configure && ${MAKE}

clean:
	${MAKE} -f psar.mk clean
	cd fsa-1.15.6 && ${MAKE} clean
