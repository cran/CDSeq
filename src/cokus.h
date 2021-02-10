#ifndef COKUS_H
#define COKUS_H

typedef unsigned long uint32;

#if (__cplusplus - 0) >= 201703L
  #define __REGISTER
#else
  #define __REGISTER register
#endif

void seedMT(uint32 seed);
uint32 randomMT();

#endif
