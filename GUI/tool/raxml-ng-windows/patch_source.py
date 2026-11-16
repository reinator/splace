"""Patches the RAxML-NG 2.0.4 source tree (argv[1]) so it builds and runs on
Windows with llvm-mingw. See build.sh."""
import sys,re,os
root=sys.argv[1]
def sub(path, old, new, count=1):
    p=os.path.join(root,path); s=open(p).read()
    if old not in s: raise SystemExit(f'not found in {path}: {old[:60]!r}')
    s=s.replace(old,new,count); open(p,'w').write(s)

# coraxlib is linked statically: no dllexport (thread-local vars can't be exported).
sub('libs/coraxlib/src/corax/core/common.h','#define CORAX_EXPORT __declspec(dllexport)','#define CORAX_EXPORT')
# Thread pinning is Linux only; std::thread handles are not pthread_t on Windows.
sub('src/ParallelContext.cpp','#ifdef _RAXML_PTHREADS\nstatic void pin_thread','#if defined(_RAXML_PTHREADS) && defined(__linux__)\nstatic void pin_thread')
sub('src/ParallelContext.cpp','''    if (opts.thread_pinning && _num_threads > 1)
    {
      auto handle = i > 0 ? _threads.back().native_handle() : pthread_self();
      pin_thread(i, handle);
    }''','''#ifdef __linux__
    if (opts.thread_pinning && _num_threads > 1)
    {
      auto handle = i > 0 ? _threads.back().native_handle() : pthread_self();
      pin_thread(i, handle);
    }
#endif''')
# unsigned long is 32-bit on Windows.
sub('src/main.cpp','std::max(2UL, total_freerate_category_count)','std::max<size_t>(2, total_freerate_category_count)')
# System utilities: Windows versions.
sub('src/util/sysutil.cpp','''#include <sys/time.h>
#include <sys/resource.h>''','''#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#include <intrin.h>
#include <io.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif''')
sub('src/util/sysutil.cpp','''  void * t = nullptr;
  int err = posix_memalign(& t, alignment, size);''','''  void * t = nullptr;
#ifdef _WIN32
  t = _aligned_malloc(size, alignment);
  int err = t ? 0 : 1;
#else
  int err = posix_memalign(& t, alignment, size);
#endif''')
sub('src/util/sysutil.cpp','#ifdef WIN32 // WINDOWS build','#ifdef _WIN32 // WINDOWS build')
sub('src/util/sysutil.cpp','''void sysutil_show_rusage()
{''','''#ifdef _WIN32
static double filetime_seconds(const FILETIME& f)
{
  ULONGLONG t = ((ULONGLONG)f.dwHighDateTime << 32) | f.dwLowDateTime;
  return t / 1.0e7;
}

void sysutil_show_rusage()
{
  FILETIME c, e, k, u;
  GetProcessTimes(GetCurrentProcess(), &c, &e, &k, &u);
  fprintf(stderr, "Time: %.3fs (user) %.3fs (sys)", filetime_seconds(u), filetime_seconds(k));
  fprintf(stderr, " Memory: %.0fMB\\n", sysutil_get_memused() * 1.0e-6);
}

unsigned long sysutil_get_memused()
{
  PROCESS_MEMORY_COUNTERS pmc;
  if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    return (unsigned long) pmc.PeakWorkingSetSize;
  return 0;
}

unsigned long sysutil_get_memtotal(bool ignore_errors)
{
  MEMORYSTATUSEX st;
  st.dwLength = sizeof(st);
  if (!GlobalMemoryStatusEx(&st))
  {
    if (ignore_errors)
      return 0;
    throw runtime_error("Cannot determine amount of RAM");
  }
  // unsigned long is 32-bit on Windows: cap at 4 GB below the limit.
  return st.ullTotalPhys > ULONG_MAX ? ULONG_MAX : (unsigned long) st.ullTotalPhys;
}
#else
void sysutil_show_rusage()
{''')
sub('src/util/sysutil.cpp','''  return si.totalram * si.mem_unit;

#endif
}''','''  return si.totalram * si.mem_unit;

#endif
}
#endif''')
sub('src/util/sysutil.cpp','''long sysutil_online_cpu_cores()
{
  return sysconf(_SC_NPROCESSORS_ONLN);''','''long sysutil_online_cpu_cores()
{
#ifdef _WIN32
  return std::thread::hardware_concurrency();
#else
  return sysconf(_SC_NPROCESSORS_ONLN);
#endif''')
sub('src/util/sysutil.cpp','''  char * real_path = realpath(path.c_str(), nullptr);''','''#ifdef _WIN32
  char * real_path = _fullpath(nullptr, path.c_str(), 0);
#else
  char * real_path = realpath(path.c_str(), nullptr);
#endif''')
# u_int32_t is a glibc name.
p=os.path.join(root,'src/util/sysutil.cpp'); t=open(p).read(); open(p,'w').write(re.sub(r'\bu_int32_t\b','uint32_t',t))
# realpath -> _fullpath into a buffer (the copy comes from our allocator).
sub('src/util/sysutil.cpp','''  char * real_path = _fullpath(nullptr, path.c_str(), 0);''','''  char full[4096];
  char * real_path = _fullpath(full, path.c_str(), sizeof(full)) ? strdup(full) : nullptr;''')
# corax_random_r returns 0..2^31-1; RAND_MAX is only 32767 on Windows, so the
# shuffle index ran past the array.
sub('libs/coraxlib/src/corax/util/stepwise.c','double r = ((double)rint / RAND_MAX);','double r = ((double)rint / 2147483648.0);')
print('patched')
