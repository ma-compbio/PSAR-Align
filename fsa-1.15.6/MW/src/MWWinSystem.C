#include "MWSystem.h"

#include <windows.h>
#include <stdio.h>
#include <direct.h>

int
MWSystem::getpid() {
  return GetCurrentProcessId();
}

static bool initialized = false;

void
MWSystem::gethostname(char *h, int len) {
  if (!initialized) {
    WSAData wsaData;
    if (WSAStartup(MAKEWORD(1, 1), &wsaData) !=0) {
      printf("Can't WSAStartup\n");
      exit(-1);
    };
    initialized = true;
  }
  
  ::gethostname(h, len);
  return;
}

double
MWSystem::gettimeofday() {
  return 0.0;
}

double
MWSystem::getcputime() {
     return 0.0;
}

void
MWSystem::sleep(int time) {
  ::Sleep(time * 1000);
}

void
MWSystem::mkdir(char *filename) {
  ::mkdir(filename);
}
