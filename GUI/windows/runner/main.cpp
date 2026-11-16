#include <flutter/dart_project.h>
#include <flutter/flutter_view_controller.h>
#include <flutter_windows.h>
#include <windows.h>

#include <algorithm>

#include "flutter_window.h"
#include "utils.h"

int APIENTRY wWinMain(_In_ HINSTANCE instance, _In_opt_ HINSTANCE prev,
                      _In_ wchar_t *command_line, _In_ int show_command) {
  // Attach to console when present (e.g., 'flutter run') or create a
  // new console when running with a debugger.
  if (!::AttachConsole(ATTACH_PARENT_PROCESS) && ::IsDebuggerPresent()) {
    CreateAndAttachConsole();
  }

  // Initialize COM, so that it is available for use in the library and/or
  // plugins.
  ::CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);

  flutter::DartProject project(L"data");

  std::vector<std::string> command_line_arguments =
      GetCommandLineArguments();

  project.set_dart_entrypoint_arguments(std::move(command_line_arguments));

  FlutterWindow window(project);
  // Open centred in the primary monitor's work area: up to 1440x900, never
  // below the 1280x720 minimum (logical pixels; the window enforces the
  // same minimum on resize).
  RECT work;
  ::SystemParametersInfo(SPI_GETWORKAREA, 0, &work, 0);
  const POINT corner = {work.left, work.top};
  const double scale =
      FlutterDesktopGetDpiForMonitor(::MonitorFromPoint(corner, MONITOR_DEFAULTTOPRIMARY)) / 96.0;
  const int work_width = static_cast<int>((work.right - work.left) / scale);
  const int work_height = static_cast<int>((work.bottom - work.top) / scale);
  const int width = std::clamp(work_width - 80, 1280, 1440);
  const int height = std::clamp(work_height - 60, 720, 900);
  Win32Window::Point origin(
      static_cast<unsigned int>(work.left / scale + std::max(0, (work_width - width) / 2)),
      static_cast<unsigned int>(work.top / scale + std::max(0, (work_height - height) / 2)));
  Win32Window::Size size(width, height);
  if (!window.Create(L"SPLACE", origin, size)) {
    return EXIT_FAILURE;
  }
  window.SetQuitOnClose(true);

  ::MSG msg;
  while (::GetMessage(&msg, nullptr, 0, 0)) {
    ::TranslateMessage(&msg);
    ::DispatchMessage(&msg);
  }

  ::CoUninitialize();
  return EXIT_SUCCESS;
}
