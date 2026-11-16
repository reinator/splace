import Cocoa
import FlutterMacOS

class MainFlutterWindow: NSWindow {
  override func awakeFromNib() {
    let flutterViewController = FlutterViewController()
    let windowFrame = self.frame
    self.contentViewController = flutterViewController
    self.setFrame(windowFrame, display: true)
    // Same minimum as the Windows and Linux runners.
    self.contentMinSize = NSSize(width: 1280, height: 720)
    if windowFrame.width < 1280 || windowFrame.height < 720 {
      self.setContentSize(NSSize(width: 1440, height: 900))
    }

    RegisterGeneratedPlugins(registry: flutterViewController)

    super.awakeFromNib()
  }
}
