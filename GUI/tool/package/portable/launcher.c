// SPLACE portable launcher (Windows).
//
// SPLACE.exe = this program + the whole app as a zip + a trailer:
//   [launcher][payload.zip][version: 32 bytes][offset: u64][size: u64]["SPLACEPK"]
// On first run the payload is extracted once (Windows' own tar.exe reads
// zip) to %LOCALAPPDATA%\SPLACE\portable\<version>, then splace.exe is
// started from there with the same arguments. Later runs start at once.
// "--splace-extract-only" extracts and exits (used to test a build).

#define WIN32_LEAN_AND_MEAN
#define UNICODE
#define _UNICODE
#include <windows.h>
#include <shlobj.h>
#include <stdio.h>
#include <stdint.h>
#include <wchar.h>

#define MAGIC "SPLACEPK"
#define TRAILER (32 + 8 + 8 + 8)

static HWND g_splash;

static LRESULT CALLBACK SplashProc(HWND h, UINT m, WPARAM w, LPARAM l) {
  if (m == WM_PAINT) {
    PAINTSTRUCT ps;
    HDC dc = BeginPaint(h, &ps);
    RECT r;
    GetClientRect(h, &r);
    HBRUSH bg = CreateSolidBrush(RGB(0x32, 0x37, 0x95));
    FillRect(dc, &r, bg);
    DeleteObject(bg);
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, RGB(255, 255, 255));
    HFONT f = CreateFontW(-18, 0, 0, 0, FW_SEMIBOLD, 0, 0, 0, DEFAULT_CHARSET, 0, 0, CLEARTYPE_QUALITY, 0, L"Segoe UI");
    HGDIOBJ old = SelectObject(dc, f);
    DrawTextW(dc, L"SPLACE\n\nPreparing the first run\u2026  /  Preparando a primeira execu\u00e7\u00e3o\u2026", -1, &r, DT_CENTER | DT_VCENTER | DT_WORDBREAK);
    SelectObject(dc, old);
    DeleteObject(f);
    EndPaint(h, &ps);
    return 0;
  }
  return DefWindowProcW(h, m, w, l);
}

static void ShowSplash(HINSTANCE inst) {
  WNDCLASSW wc = {0};
  wc.lpfnWndProc = SplashProc;
  wc.hInstance = inst;
  wc.lpszClassName = L"SplaceSplash";
  wc.hCursor = LoadCursor(NULL, IDC_WAIT);
  wc.hIcon = LoadIconW(inst, MAKEINTRESOURCEW(1));
  RegisterClassW(&wc);
  int w = 520, h = 150;
  int x = (GetSystemMetrics(SM_CXSCREEN) - w) / 2, y = (GetSystemMetrics(SM_CYSCREEN) - h) / 2;
  g_splash = CreateWindowExW(WS_EX_TOPMOST | WS_EX_TOOLWINDOW, L"SplaceSplash", L"SPLACE", WS_POPUP | WS_VISIBLE, x, y, w, h, NULL, NULL, inst, NULL);
  UpdateWindow(g_splash);
}

// Waits for [proc] while keeping the splash window responsive.
static DWORD WaitPumping(HANDLE proc) {
  for (;;) {
    DWORD r = MsgWaitForMultipleObjects(1, &proc, FALSE, INFINITE, QS_ALLINPUT);
    if (r == WAIT_OBJECT_0) break;
    MSG msg;
    while (PeekMessageW(&msg, NULL, 0, 0, PM_REMOVE)) {
      TranslateMessage(&msg);
      DispatchMessageW(&msg);
    }
  }
  DWORD code = 1;
  GetExitCodeProcess(proc, &code);
  return code;
}

static int Fail(const wchar_t *what) {
  wchar_t msg[1024];
  swprintf(msg, 1024, L"SPLACE could not start:\n%ls\n\n(error %lu)", what, GetLastError());
  if (g_splash) DestroyWindow(g_splash);
  MessageBoxW(NULL, msg, L"SPLACE", MB_ICONERROR);
  return 1;
}

static int FileExists(const wchar_t *p) { return GetFileAttributesW(p) != INVALID_FILE_ATTRIBUTES; }

int WINAPI wWinMain(HINSTANCE inst, HINSTANCE prev, PWSTR cmdLine, int show) {
  (void)prev;
  (void)show;
  wchar_t self[MAX_PATH * 4];
  if (!GetModuleFileNameW(NULL, self, MAX_PATH * 4)) return Fail(L"own path");

  FILE *f = _wfopen(self, L"rb");
  if (!f) return Fail(L"cannot read the program file");
  _fseeki64(f, -TRAILER, SEEK_END);
  char version[33] = {0};
  uint64_t offset = 0, size = 0;
  char magic[9] = {0};
  if (fread(version, 1, 32, f) != 32 || fread(&offset, 8, 1, f) != 1 || fread(&size, 8, 1, f) != 1 || fread(magic, 1, 8, f) != 8 ||
      memcmp(magic, MAGIC, 8) != 0) {
    fclose(f);
    return Fail(L"this file has no SPLACE package inside (incomplete download?)");
  }

  // %LOCALAPPDATA%\SPLACE\portable\<version>
  wchar_t base[MAX_PATH], dir[MAX_PATH * 2], app[MAX_PATH * 2], ok[MAX_PATH * 2], zip[MAX_PATH * 2];
  if (FAILED(SHGetFolderPathW(NULL, CSIDL_LOCAL_APPDATA, NULL, 0, base))) return Fail(L"LOCALAPPDATA");
  wchar_t wver[33];
  MultiByteToWideChar(CP_UTF8, 0, version, -1, wver, 33);
  swprintf(dir, MAX_PATH * 2, L"%ls\\SPLACE\\portable\\%ls", base, wver);
  swprintf(app, MAX_PATH * 2, L"%ls\\splace.exe", dir);
  swprintf(ok, MAX_PATH * 2, L"%ls\\.extracted", dir);
  swprintf(zip, MAX_PATH * 2, L"%ls\\payload.zip", dir);
  int extractOnly = cmdLine && wcsstr(cmdLine, L"--splace-extract-only") != NULL;

  if (!FileExists(ok)) {
    if (!extractOnly) ShowSplash(inst);
    SHCreateDirectoryExW(NULL, dir, NULL);
    // The payload to a zip next to the app, then tar.exe extracts it.
    FILE *out = _wfopen(zip, L"wb");
    if (!out) {
      fclose(f);
      return Fail(dir);
    }
    _fseeki64(f, (int64_t)offset, SEEK_SET);
    static char buf[1 << 20];
    uint64_t left = size;
    while (left > 0) {
      size_t n = (size_t)(left < sizeof buf ? left : sizeof buf);
      if (fread(buf, 1, n, f) != n || fwrite(buf, 1, n, out) != n) {
        fclose(out);
        fclose(f);
        return Fail(L"could not copy the package");
      }
      left -= n;
    }
    fclose(out);
    fclose(f);

    wchar_t tar[MAX_PATH], cmd[MAX_PATH * 6];
    GetSystemDirectoryW(tar, MAX_PATH);
    swprintf(cmd, MAX_PATH * 6, L"\"%ls\\tar.exe\" -xf \"%ls\" -C \"%ls\"", tar, zip, dir);
    STARTUPINFOW si = {sizeof si};
    PROCESS_INFORMATION pi;
    if (!CreateProcessW(NULL, cmd, NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, dir, &si, &pi)) return Fail(L"tar.exe (Windows 10 1803 or newer is needed)");
    DWORD code = WaitPumping(pi.hProcess);
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    DeleteFileW(zip);
    if (code != 0 || !FileExists(app)) return Fail(L"extraction failed");
    HANDLE h = CreateFileW(ok, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_HIDDEN, NULL);
    if (h != INVALID_HANDLE_VALUE) CloseHandle(h);
    // Older versions unpacked before are removed (best effort).
    wchar_t pattern[MAX_PATH * 2], parent[MAX_PATH * 2];
    swprintf(parent, MAX_PATH * 2, L"%ls\\SPLACE\\portable", base);
    swprintf(pattern, MAX_PATH * 2, L"%ls\\*", parent);
    WIN32_FIND_DATAW fd;
    HANDLE find = FindFirstFileW(pattern, &fd);
    if (find != INVALID_HANDLE_VALUE) {
      do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) || fd.cFileName[0] == L'.' || wcscmp(fd.cFileName, wver) == 0) continue;
        wchar_t rm[MAX_PATH * 6];
        swprintf(rm, MAX_PATH * 6, L"cmd.exe /c rmdir /s /q \"%ls\\%ls\"", parent, fd.cFileName);
        STARTUPINFOW rsi = {sizeof rsi};
        PROCESS_INFORMATION rpi;
        if (CreateProcessW(NULL, rm, NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, NULL, &rsi, &rpi)) {
          CloseHandle(rpi.hThread);
          CloseHandle(rpi.hProcess);
        }
      } while (FindNextFileW(find, &fd));
      FindClose(find);
    }
  } else {
    fclose(f);
  }
  if (g_splash) DestroyWindow(g_splash);
  if (extractOnly) return 0;

  // splace.exe with the same arguments, from its own folder.
  wchar_t *cmd = (wchar_t *)HeapAlloc(GetProcessHeap(), 0, sizeof(wchar_t) * (wcslen(app) + wcslen(cmdLine ? cmdLine : L"") + 8));
  swprintf(cmd, wcslen(app) + wcslen(cmdLine ? cmdLine : L"") + 8, L"\"%ls\" %ls", app, cmdLine ? cmdLine : L"");
  STARTUPINFOW si = {sizeof si};
  PROCESS_INFORMATION pi;
  if (!CreateProcessW(NULL, cmd, NULL, NULL, FALSE, 0, NULL, dir, &si, &pi)) return Fail(app);
  CloseHandle(pi.hThread);
  CloseHandle(pi.hProcess);
  return 0;
}
