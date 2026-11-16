/* Helpers missing from the Windows C runtime, for the SPLACE build of
   RAxML-NG (force-included in every file). */
#ifndef SPLACE_WIN_COMPAT_H
#define SPLACE_WIN_COMPAT_H
#ifdef _WIN32
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#ifdef __cplusplus
extern "C" {
#endif
static inline int splace_vasprintf(char **out, const char *fmt, va_list ap)
{
  va_list cp;
  va_copy(cp, ap);
  int n = vsnprintf(NULL, 0, fmt, cp);
  va_end(cp);
  if (n < 0) { *out = NULL; return -1; }
  *out = (char *) malloc((size_t) n + 1);
  if (!*out) return -1;
  return vsnprintf(*out, (size_t) n + 1, fmt, ap);
}
static inline int splace_asprintf(char **out, const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  int n = splace_vasprintf(out, fmt, ap);
  va_end(ap);
  return n;
}
#ifdef __cplusplus
}
#endif
/* glibc's rand()/srand() (TYPE_3 additive feedback generator, values up to
   2^31-1), so a seed gives the same analysis as on Linux; the Windows
   rand() only reaches 32767. */
static int splace_rand_i = -1;
static long long splace_rand_tbl[31];
static int splace_rand_f, splace_rand_b;
static inline void splace_srand(unsigned int seed)
{
  int i;
  if (seed == 0) seed = 1;
  splace_rand_tbl[0] = (int) seed;
  for (i = 1; i < 31; i++) {
    long long word = (16807LL * splace_rand_tbl[i - 1]) % 2147483647LL;
    if (word < 0) word += 2147483647LL;
    splace_rand_tbl[i] = word;
  }
  splace_rand_f = 3; splace_rand_b = 0;
  splace_rand_i = 0;
  /* discard the first 310 outputs, as glibc does */
  for (i = 0; i < 310; i++) {
    unsigned int v = (unsigned int) splace_rand_tbl[splace_rand_f] + (unsigned int) splace_rand_tbl[splace_rand_b];
    splace_rand_tbl[splace_rand_f] = (int) v;
    if (++splace_rand_f >= 31) splace_rand_f = 0;
    if (++splace_rand_b >= 31) splace_rand_b = 0;
  }
}
static inline int splace_rand(void)
{
  if (splace_rand_i < 0) splace_srand(1);
  unsigned int v = (unsigned int) splace_rand_tbl[splace_rand_f] + (unsigned int) splace_rand_tbl[splace_rand_b];
  splace_rand_tbl[splace_rand_f] = (int) v;
  int result = (int) (v >> 1);
  if (++splace_rand_f >= 31) splace_rand_f = 0;
  if (++splace_rand_b >= 31) splace_rand_b = 0;
  return result;
}
#ifdef __cplusplus
extern "C" {
#endif
void *splace_aligned_malloc(size_t n, size_t align);
char *splace_strdup(const char *s);
#ifdef __cplusplus
}
#endif
#define _aligned_malloc(n, a) splace_aligned_malloc((n), (a))
#define _aligned_free(p) free(p)
#define strdup splace_strdup
#define _strdup splace_strdup
#define asprintf splace_asprintf
#define rand splace_rand
#define srand splace_srand
#define vasprintf splace_vasprintf
#endif
#endif
