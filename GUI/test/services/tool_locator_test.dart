import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/services/tool_locator.dart';

void main() {
  test('Windows paths become WSL paths, also after key=', () {
    expect(windowsToWslPath(r'C:\Users\x\run\m.phy'), '/mnt/c/Users/x/run/m.phy');
    expect(windowsToWslPath(r'auto=D:\p\partitions_raxml.txt'), 'auto=/mnt/d/p/partitions_raxml.txt');
    expect(windowsToWslPath('criterion=bic'), 'criterion=bic');
    expect(windowsToWslPath('GTR+G'), 'GTR+G');
  });
}
