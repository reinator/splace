# Turns ClipKIT, PhyKIT and GMYC into standalone programs (PyInstaller, one
# folder each: <Dest>\<name>\<name>.exe), so the users need no Python. One
# folder each: they start at once, a one-file build unpacks itself on every
# call. The app finds them in its tools\ folder before any Python.
#
#   powershell -ExecutionPolicy Bypass -File tool\package\python_tools.ps1            # into GUI\tools\win (every build copies it)
#   powershell -ExecutionPolicy Bypass -File tool\package\python_tools.ps1 -Dest <dir>
#
# Needs Python 3. The virtual environment is kept in build\pytools.
param([string]$Dest = '')
$ErrorActionPreference = 'Stop'
$root = Resolve-Path "$PSScriptRoot\..\.."
if (-not $Dest) { $Dest = "$root\tools\win" }
New-Item -ItemType Directory -Force $Dest | Out-Null

$venv = "$root\build\pytools"
if (-not (Test-Path "$venv\Scripts\python.exe")) {
  $py = (Get-Command python -ErrorAction SilentlyContinue).Source
  if (-not $py) { throw 'Python 3 not found (needed to build ClipKIT, PhyKIT and GMYC)' }
  & $py -m venv $venv
}
$python = "$venv\Scripts\python.exe"
& $python -m pip install --upgrade pip pyinstaller clipkit phykit ete3 scipy six
if ($LASTEXITCODE -ne 0) { throw 'pip install failed' }

foreach ($t in @(@{ Name = 'clipkit'; Module = 'clipkit.clipkit' }, @{ Name = 'phykit'; Module = 'phykit.phykit' })) {
  $entry = "$venv\$($t.Name)_entry.py"
  "import sys`nfrom $($t.Module) import main`nsys.exit(main())" | Set-Content $entry -Encoding utf8
  & $python -m PyInstaller --noconfirm --onedir --console --name $t.Name --collect-all $t.Name --collect-data Bio `
    --distpath $Dest --workpath "$venv\work" --specpath $venv $entry
  if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed for $($t.Name)" }
}
# GMYC (pGMYC, the iTaxoTools Python 3 port in tool\gmyc, with SPLACE's
# command line): ete3 and scipy only, no Qt or matplotlib.
& $python -m PyInstaller --noconfirm --onedir --console --name gmyc --paths "$root\tool\gmyc" --hidden-import GMYC `
  --collect-submodules ete3 --exclude-module PyQt5 --exclude-module matplotlib `
  --distpath $Dest --workpath "$venv\work" --specpath $venv "$root\tool\gmyc\splace_gmyc.py"
if ($LASTEXITCODE -ne 0) { throw 'PyInstaller failed for GMYC' }

# Each one must start on its own (ClipKIT exits with 1 after --version, so
# the output is checked, not the exit code).
foreach ($n in 'clipkit', 'phykit', 'gmyc') {
  $exe = "$Dest\$n\$n.exe"
  $v = (& $exe --version 2>&1 | Out-String).Trim()
  if (-not $v -or $v -match 'Traceback|ModuleNotFoundError') { throw "$n does not start: $v" }
  Write-Host "$n -> $exe ($(($v -split "`n")[0]))"
}
