/* SPLACE build of RAxML-NG for Windows: one allocator for the whole
   program. coraxlib allocates aligned memory (_aligned_malloc on Windows)
   and sometimes releases it with free(), which glibc accepts and the
   Windows CRT does not. Here every block is 64-byte aligned and free()
   takes any of them. */
#include <windows.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>

#define SPLACE_ALIGN 64
#define SPLACE_HDR (2 * sizeof(void *))

static void *splace_block(size_t n, size_t align)
{
  if (n == 0) n = 1;
  if (align < SPLACE_ALIGN) align = SPLACE_ALIGN;
  if (n > SIZE_MAX - align - SPLACE_HDR) { errno = ENOMEM; return NULL; }
  char *raw = (char *) HeapAlloc(GetProcessHeap(), 0, n + align + SPLACE_HDR);
  if (!raw) { errno = ENOMEM; return NULL; }
  uintptr_t a = ((uintptr_t) raw + SPLACE_HDR + align - 1) & ~(uintptr_t) (align - 1);
  ((void **) a)[-1] = raw;
  ((size_t *) a)[-2] = n;
  return (void *) a;
}

void *malloc(size_t n) { return splace_block(n, SPLACE_ALIGN); }

void free(void *p)
{
  if (p) HeapFree(GetProcessHeap(), 0, ((void **) p)[-1]);
}

void *calloc(size_t c, size_t s)
{
  if (s && c > SIZE_MAX / s) { errno = ENOMEM; return NULL; }
  void *p = malloc(c * s);
  if (p) memset(p, 0, c * s);
  return p;
}

void *realloc(void *p, size_t n)
{
  if (!p) return malloc(n);
  if (n == 0) { free(p); return NULL; }
  size_t old = ((size_t *) p)[-2];
  void *q = malloc(n);
  if (!q) return NULL;
  memcpy(q, p, old < n ? old : n);
  free(p);
  return q;
}

void *splace_aligned_malloc(size_t n, size_t align) { return splace_block(n, align); }

char *splace_strdup(const char *s)
{
  size_t n = strlen(s) + 1;
  char *d = (char *) malloc(n);
  if (d) memcpy(d, s, n);
  return d;
}
