/// Global pairwise alignment of two DNA stretches (Needleman–Wunsch, linear
/// gap penalty): used to draw the alignment of a BOLD hit, whose aligned
/// sequences the ID engine does not return (only where they align).
library;

const _match = 2, _mismatch = -3, _gap = -5;

/// [a] and [b] aligned end to end, with '-' in the gaps.
(String, String) alignPair(String a, String b) {
  a = a.toUpperCase();
  b = b.toUpperCase();
  final n = a.length, m = b.length;
  final w = m + 1;
  final score = List<int>.filled((n + 1) * w, 0);
  // 0 diagonal, 1 up (gap in b), 2 left (gap in a).
  final move = List<int>.filled((n + 1) * w, 0);
  for (var i = 1; i <= n; i++) {
    score[i * w] = i * _gap;
    move[i * w] = 1;
  }
  for (var j = 1; j <= m; j++) {
    score[j] = j * _gap;
    move[j] = 2;
  }
  for (var i = 1; i <= n; i++) {
    final ai = a.codeUnitAt(i - 1);
    for (var j = 1; j <= m; j++) {
      final bj = b.codeUnitAt(j - 1);
      // N matches anything without reward or penalty.
      final s = ai == 78 || bj == 78 ? 0 : (ai == bj ? _match : _mismatch);
      final d = score[(i - 1) * w + j - 1] + s;
      final u = score[(i - 1) * w + j] + _gap;
      final l = score[i * w + j - 1] + _gap;
      var best = d, mv = 0;
      if (u > best) {
        best = u;
        mv = 1;
      }
      if (l > best) {
        best = l;
        mv = 2;
      }
      score[i * w + j] = best;
      move[i * w + j] = mv;
    }
  }
  final x = StringBuffer(), y = StringBuffer();
  var i = n, j = m;
  final ra = <String>[], rb = <String>[];
  while (i > 0 || j > 0) {
    final mv = move[i * w + j];
    if (i > 0 && j > 0 && mv == 0) {
      ra.add(a[i - 1]);
      rb.add(b[j - 1]);
      i--;
      j--;
    } else if (i > 0 && (j == 0 || mv == 1)) {
      ra.add(a[i - 1]);
      rb.add('-');
      i--;
    } else {
      ra.add('-');
      rb.add(b[j - 1]);
      j--;
    }
  }
  for (var k = ra.length - 1; k >= 0; k--) {
    x.write(ra[k]);
    y.write(rb[k]);
  }
  return (x.toString(), y.toString());
}
