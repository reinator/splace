/* Windows helpers for the SPLACE builds of ASAP and ABGD (force-included in
   every file with -include). Elsewhere it does nothing. */
#ifndef SPLACE_DELIM_COMPAT_H
#define SPLACE_DELIM_COMPAT_H
#ifdef _WIN32
#include <stdlib.h>
#include <string.h>
#include <io.h>
#include <direct.h>
#include <getopt.h>
#include <process.h>
/* POSIX mkdir(path, mode): the Windows one takes only the path. */
#define mkdir(path, mode) _mkdir(path)
#define bzero(p, n) memset((p), 0, (n))
/* glibc's rand()/srand() (values up to 2^31-1, as on Linux and macOS): the
   Windows rand() stops at 32767, which coarsens ASAP's permutation tests. */
int splace_rand(void);
void splace_srand(unsigned int seed);
#undef RAND_MAX
#define RAND_MAX 2147483647
#define rand() splace_rand()
#define srand(s) splace_srand(s)
#endif
#endif
