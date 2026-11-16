# RAxML-NG for Windows

RAxML-NG has no official Windows release. SPLACE ships `raxml-ng.exe`
built from the official RAxML-NG 2.0.4 source with the patches in this
folder (`build.sh` reproduces it on Linux or WSL). RAxML-NG is licensed under
the GNU AGPL v3: the source is the official 2.0.4 release, and the
changes are `patch_source.py`, `win_compat.h` and `splace_alloc.c`.

Checked against the official Linux binary: the same models chosen by MOOSE,
the same ML topology and log-likelihoods equal to the 4th decimal.
Bootstrap replicates differ from Linux for the same seed, because the C++
library draws its random numbers in a different way.
