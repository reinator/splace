import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/services/ncbi.dart';

/// The live GBIF backbone (SPLACE_NETWORK=1).
void main() {
  final network = Platform.environment['SPLACE_NETWORK'] == '1';
  test('synonyms resolve to the accepted species; unknown names to null', () async {
    final ncbi = NcbiClient(cacheDir: Directory.systemTemp.path);
    expect(await ncbi.gbifAccepted('Bufo marinus'), 'Rhinella marina');
    expect(await ncbi.gbifAccepted('Chaunus marinus'), 'Rhinella marina');
    expect(await ncbi.gbifAccepted('Rhinella marina'), 'Rhinella marina');
    expect(await ncbi.gbifAccepted('Qwertyus asdfgh'), isNull);
  }, skip: network ? false : 'set SPLACE_NETWORK=1');
}
