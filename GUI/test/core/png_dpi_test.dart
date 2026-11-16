import 'dart:typed_data';

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/ui/widgets/tree_viewer.dart';

void main() {
  test('pngWithDpi inserts a pHYs chunk after IHDR', () {
    // Signature + a 13-byte IHDR chunk + IEND (content does not matter).
    final png = Uint8List.fromList([
      137, 80, 78, 71, 13, 10, 26, 10, //
      0, 0, 0, 13, 73, 72, 68, 82, ...List.filled(13, 0), 0, 0, 0, 0, //
      0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130,
    ]);
    final out = pngWithDpi(png, 300);
    expect(out.length, png.length + 21);
    final b = ByteData.sublistView(out);
    expect(b.getUint32(33), 9);
    expect(String.fromCharCodes(out.sublist(37, 41)), 'pHYs');
    expect(b.getUint32(41), 11811); // 300 dpi in pixels per metre
    expect(b.getUint32(45), 11811);
    expect(out[49], 1);
    expect(b.getUint32(50), 0x78A53F76);
    expect(out.sublist(54), png.sublist(33));
  });
}
