AM_CPPFLAGS = -I$(top_srcdir)/src -Wno-deprecated

noinst_LIBRARIES = libmanager.a

libmanager_a_SOURCES = \
	manager.cc \
	mw_adapter.cc \
	db_adapter.cc \
	db_misc.cc 

noinst_HEADERS = \
	manager.h \
	mw_adapter.h \
	db_adapter.h \
	db_misc.h

if HAVE_CONDOR

AM_CPPFLAGS += -I$(top_srcdir)/MW/src -I$(top_srcdir)/MW/src/MWControlTasks -I$(top_srcdir)/MW/src/RMComm -I$(top_srcdir)/MW/src/RMComm/MW-Socket

libmanager_a_SOURCES += \
	mw_master.cc \
	mw_worker.cc \
	mw_task.cc \
	transfer_data.cc

noinst_HEADERS += \
	mw_master.h \
	mw_worker.h \
	mw_task.h

else

if HAVE_POSTGRES
libmanager_a_SOURCES += \
	transfer_data.cc 
endif

endif

if HAVE_POSTGRES

libmanager_a_SOURCES += \
	db_postgres.cc 

noinst_HEADERS += \
	db_postgres.h 

endif
