if HAVE_CONDOR
  MAYBE_MW = MW
endif

ALWAYS_BUILT = src/util src/math src/seq src/manager src/fsa src/annealing src/main perl examples doc html tests
SUBDIRS = $(MAYBE_MW) $(ALWAYS_BUILT)
DIST_SUBDIRS = $(ALWAYS_BUILT)

EXTRA_DIST = LICENSE COPYING $(MW_DIST)

# Possibly regenerate version file if this is a git repository.
$(top_srcdir)/version.m4: FORCE
	@cd $(top_srcdir); if test -r version.sh; then ./version.sh; fi
.PHONY: FORCE

MW_DIST = \
	MW/acinclude.m4 \
	MW/aclocal.m4 \
	MW/configure \
	MW/configure.in \
	MW/COPYING \
	MW/doxygen.conf \
	MW/INSTALL \
	MW/install-sh \
	MW/Makefile.in \
	MW/missing \
	MW/mkinstalldirs \
	MW/README \
	MW/src/BlackBox/MWDriver_blackbox.C \
	MW/src/BlackBox/MWDriver_blackbox.h \
	MW/src/BlackBox/MWTask_blackbox.C \
	MW/src/BlackBox/MWTask_blackbox.h \
	MW/src/BlackBox/MWWorkerMain_blackbox.C \
	MW/src/BlackBox/MWWorker_blackbox.C \
	MW/src/BlackBox/MWWorker_blackbox.h \
	MW/src/BlackBox/Makefile.in \
	MW/src/MW.C \
	MW/src/MW.h \
	MW/src/MWControlTasks/MWNWSTask.C \
	MW/src/MWControlTasks/MWNWSTask.h \
	MW/src/MWControlTasks/Makefile.in \
	MW/src/MWDriver.C \
	MW/src/MWDriver.h \
	MW/src/MWGroup.C \
	MW/src/MWGroup.h \
	MW/src/MWList.C \
	MW/src/MWList.h \
	MW/src/MWMasterMain.C \
	MW/src/MWStats.C \
	MW/src/MWStats.h \
	MW/src/MWSystem.h \
	MW/src/MWTask.C \
	MW/src/MWTask.h \
	MW/src/MWTaskContainer.C \
	MW/src/MWTaskContainer.h \
	MW/src/MWUnixSystem.C \
	MW/src/MWWinSystem.C \
	MW/src/MWWorker.C \
	MW/src/MWWorker.h \
	MW/src/MWWorkerID.C \
	MW/src/MWWorkerID.h \
	MW/src/MWWorkerMain.C \
	MW/src/MWprintf.C \
	MW/src/MWprintf.h \
	MW/src/Makefile.in \
	MW/src/RMComm/MW-CondorPVM/MWCondorPvmRC.C \
	MW/src/RMComm/MW-CondorPVM/MWCondorPvmRC.h \
	MW/src/RMComm/MW-CondorPVM/Makefile.in \
	MW/src/RMComm/MW-File/MWFileError.h \
	MW/src/RMComm/MW-File/MWFileRC.C \
	MW/src/RMComm/MW-File/MWFileRC.h \
	MW/src/RMComm/MW-File/MWFileRCSymbol.h \
	MW/src/RMComm/MW-File/MWFileSend.h \
	MW/src/RMComm/MW-File/MWFileTypes.h \
	MW/src/RMComm/MW-File/MWFileWorker.h \
	MW/src/RMComm/MW-File/Makefile.in \
	MW/src/RMComm/MW-File/chirp_client.c \
	MW/src/RMComm/MW-File/chirp_client.h \
	MW/src/RMComm/MW-File/chirp_protocol.h \
	MW/src/RMComm/MW-Independent/MWIndRC.C \
	MW/src/RMComm/MW-Independent/MWIndRC.h \
	MW/src/RMComm/MW-Independent/Makefile.in \
	MW/src/RMComm/MW-Socket/MWSocketRC.C \
	MW/src/RMComm/MW-Socket/MWSocketRC.h \
	MW/src/RMComm/MW-Socket/Makefile.in \
	MW/src/RMComm/MW-StaticMPI/MWStaticMPIRC.C \
	MW/src/RMComm/MW-StaticMPI/MWStaticMPIRC.h \
	MW/src/RMComm/MW-StaticMPI/Makefile.in \
	MW/src/RMComm/MWRMComm.C \
	MW/src/RMComm/MWRMComm.h \
	MW/src/RMComm/Makefile.in
