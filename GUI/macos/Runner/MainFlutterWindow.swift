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

    // The title bar shows the project file (with its proxy icon) and
    // whether it has unsaved changes (the dot in the close button).
    let channel = FlutterMethodChannel(name: "splace/window", binaryMessenger: flutterViewController.engine.binaryMessenger)
    channel.setMethodCallHandler { [weak self] call, result in
      guard call.method == "setTitle", let args = call.arguments as? [String: Any] else {
        result(FlutterMethodNotImplemented)
        return
      }
      if let title = args["title"] as? String { self?.title = title }
      if let file = args["file"] as? String {
        self?.representedURL = URL(fileURLWithPath: file)
      } else {
        self?.representedURL = nil
      }
      self?.isDocumentEdited = (args["edited"] as? Bool) ?? false
      result(nil)
    }

    super.awakeFromNib()
  }
}
