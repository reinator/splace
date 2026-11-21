/* glibc's TYPE_3 additive feedback generator (rand()/srand()), for the
   Windows builds of ASAP and ABGD; see splace_compat.h. */
static long long tbl[31];
static int f_i, b_i, ready;

void splace_srand(unsigned int seed)
{
  int i;
  if (seed == 0) seed = 1;
  tbl[0] = (int) seed;
  for (i = 1; i < 31; i++) {
    long long word = (16807LL * tbl[i - 1]) % 2147483647LL;
    if (word < 0) word += 2147483647LL;
    tbl[i] = word;
  }
  f_i = 3;
  b_i = 0;
  ready = 1;
  for (i = 0; i < 310; i++) {
    unsigned int v = (unsigned int) tbl[f_i] + (unsigned int) tbl[b_i];
    tbl[f_i] = (int) v;
    if (++f_i >= 31) f_i = 0;
    if (++b_i >= 31) b_i = 0;
  }
}

int splace_rand(void)
{
  unsigned int v;
  int r;
  if (!ready) splace_srand(1);
  v = (unsigned int) tbl[f_i] + (unsigned int) tbl[b_i];
  tbl[f_i] = (int) v;
  r = (int) (v >> 1);
  if (++f_i >= 31) f_i = 0;
  if (++b_i >= 31) b_i = 0;
  return r;
}
