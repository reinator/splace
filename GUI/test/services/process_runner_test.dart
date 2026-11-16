import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/services/process_runner.dart';

void main() {
  test('crash exit codes are explained', () {
    // RAxML-NG on Windows: 0xC0000005, reported signed or unsigned.
    expect(exitCodeMeaning(-1073741819), contains('0xC0000005: the program crashed (memory access violation)'));
    expect(exitCodeMeaning(3221225477), contains('memory access violation'));
    expect(exitCodeMeaning(-1073741571), contains('stack overflow'));
    expect(exitCodeMeaning(137), contains('out of memory'));
    expect(exitCodeMeaning(2), 'exit code 2');
  });
}
