# Builds SPLACE for Windows (x64) as ONE portable program:
#   dist\SPLACE-<version>-windows-x64.exe
# No installation: copy it anywhere and double-click. It carries the app,
# the Flutter runtime, MAFFT, trimAl, IQ-TREE, RAxML-NG, ClipKIT, PhyKIT
# and the Visual C++ runtime; the first run unpacks them once to
# %LOCALAPPDATA%\SPLACE\portable\<version> (later runs start at once).
#
#   powershell -ExecutionPolicy Bypass -File tool\package\windows.ps1
#
# Needs Flutter in PATH, the Visual Studio C++ build tools (already needed
# by Flutter for Windows) and Python 3 (ClipKIT and PhyKIT are turned into
# standalone programs with PyInstaller, so the users need no Python).
#   -SkipBuild       pack the last release build again
#   -NoPythonTools   leave ClipKIT and PhyKIT out (the app installs them with pip)
#   -Installer       also make an Inno Setup installer (needs Inno Setup 6)
param([switch]$SkipBuild, [switch]$NoPythonTools, [switch]$Installer)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path "$PSScriptRoot\..\..").Path
Set-Location $root
$version = (Select-String -Path 'lib\src\ui\common.dart' -Pattern "appVersion = '([^']+)'").Matches[0].Groups[1].Value

if (-not $SkipBuild) {
  flutter build windows --release
  if ($LASTEXITCODE -ne 0) { throw 'flutter build windows failed' }
}
$bundle = "$root\build\windows\x64\runner\Release"
if (-not (Test-Path "$bundle\splace.exe")) { throw "No release build in $bundle" }
if (-not (Test-Path "$bundle\tools\raxml-ng.exe")) { Write-Warning "tools\ is missing from $bundle (electron\bin\win not found?)" }

# ClipKIT and PhyKIT as standalone programs in tools\<name>\<name>.exe
# (one folder each: they start at once, a one-file build unpacks itself on
# every call). The app finds them there before any Python.
if (-not $NoPythonTools) {
  $py = (Get-Command python -ErrorAction SilentlyContinue).Source
  if (-not $py) { throw 'Python 3 not found (needed for ClipKIT and PhyKIT; -NoPythonTools skips them)' }
  $venv = "$root\build\pytools"
  if (-not (Test-Path "$venv\Scripts\python.exe")) { & $py -m venv $venv }
  & "$venv\Scripts\python.exe" -m pip install --upgrade pip pyinstaller clipkit phykit
  if ($LASTEXITCODE -ne 0) { throw 'pip install failed' }
  foreach ($t in @(@{ Name = 'clipkit'; Module = 'clipkit.clipkit' }, @{ Name = 'phykit'; Module = 'phykit.phykit' })) {
    $entry = "$venv\$($t.Name)_entry.py"
    "import sys`nfrom $($t.Module) import main`nsys.exit(main())" | Set-Content $entry -Encoding utf8
    & "$venv\Scripts\pyinstaller.exe" --noconfirm --onedir --console --name $t.Name --collect-all $t.Name --collect-data Bio `
      --distpath "$bundle\tools" --workpath "$venv\work" --specpath $venv $entry
    if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed for $($t.Name)" }
  }
}

# The Visual C++ runtime goes next to splace.exe, for computers without the
# VC++ redistributable installed.
foreach ($dll in 'msvcp140.dll', 'vcruntime140.dll', 'vcruntime140_1.dll') {
  $src = Join-Path $env:WINDIR "System32\$dll"
  if (Test-Path $src) { Copy-Item $src $bundle -Force } else { Write-Warning "$dll not found" }
}

New-Item -ItemType Directory -Force "$root\dist" | Out-Null
$work = "$root\build\portable"
New-Item -ItemType Directory -Force $work | Out-Null

# 1. The launcher (C, static runtime, with the app icon).
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vs = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vs) { throw 'Visual Studio C++ build tools not found' }
$vcvars = "$vs\VC\Auxiliary\Build\vcvars64.bat"
$icon = "$root\windows\runner\resources\app_icon.ico"
"1 ICON `"$($icon -replace '\\', '\\')`"" | Set-Content "$work\launcher.rc" -Encoding ascii
$src = "$PSScriptRoot\portable\launcher.c"
cmd /c "`"$vcvars`" >nul && cd /d `"$work`" && rc /nologo launcher.rc && cl /nologo /O2 /MT /W3 `"$src`" launcher.res /Fe:launcher.exe /link /SUBSYSTEM:WINDOWS user32.lib gdi32.lib shell32.lib"
if ($LASTEXITCODE -ne 0) { throw 'Could not compile the launcher' }

# 2. The app as a zip.
$zip = "$work\payload.zip"
Remove-Item $zip -ErrorAction SilentlyContinue
& "$env:WINDIR\System32\tar.exe" -a -c -f $zip -C $bundle *
if ($LASTEXITCODE -ne 0) { throw 'Could not zip the app' }

# 3. launcher + zip + trailer (version, offset, size, "SPLACEPK").
$tag = "$version-$((Get-Date).ToString('yyyyMMddHHmm'))"
$exe = "$root\dist\SPLACE-$version-windows-x64.exe"
$launcher = [IO.File]::ReadAllBytes("$work\launcher.exe")
$out = [IO.File]::Create($exe)
try {
  $out.Write($launcher, 0, $launcher.Length)
  $offset = [uint64]$out.Position
  $in = [IO.File]::OpenRead($zip)
  try { $in.CopyTo($out) } finally { $in.Dispose() }
  $size = [uint64]$out.Position - $offset
  $ver = New-Object byte[] 32
  $vb = [Text.Encoding]::UTF8.GetBytes($tag)
  [Array]::Copy($vb, $ver, [Math]::Min($vb.Length, 31))
  $out.Write($ver, 0, 32)
  $out.Write([BitConverter]::GetBytes($offset), 0, 8)
  $out.Write([BitConverter]::GetBytes($size), 0, 8)
  $magic = [Text.Encoding]::ASCII.GetBytes('SPLACEPK')
  $out.Write($magic, 0, 8)
} finally { $out.Dispose() }
Write-Host ("Portable program: $exe ({0:N0} MB)" -f ((Get-Item $exe).Length / 1MB))

if ($Installer) {
  $iscc = (Get-Command iscc.exe -ErrorAction SilentlyContinue).Source
  foreach ($c in "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe", "$env:ProgramFiles\Inno Setup 6\ISCC.exe", "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe") {
    if (-not $iscc -and (Test-Path $c)) { $iscc = $c }
  }
  if (-not $iscc) { throw 'Inno Setup 6 not found. Install it with:  winget install --id JRSoftware.InnoSetup -e' }
  & $iscc /Q "/DAppVersion=$version" "/DBundle=$bundle" "/DOutDir=$root\dist" "/DIcon=$icon" "$PSScriptRoot\splace.iss"
  if ($LASTEXITCODE -ne 0) { throw 'Inno Setup failed' }
  Write-Host "Installer: $root\dist\SPLACE-$version-windows-x64-setup.exe"
}
