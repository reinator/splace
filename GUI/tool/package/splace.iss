; Inno Setup script for SPLACE (called by windows.ps1 with /DAppVersion,
; /DBundle, /DOutDir and /DIcon). Installs for the current user only, so no
; administrator rights are needed and the tools can write next to the app.

[Setup]
AppId={{6C1E6A43-3B0F-4C8B-9E55-2C7C2A3F5D10}
AppName=SPLACE
AppVersion={#AppVersion}
AppVerName=SPLACE {#AppVersion}
AppPublisher=ITV Genomics
AppPublisherURL=https://github.com/luanrabelo/splace
DefaultDirName={localappdata}\Programs\SPLACE
DefaultGroupName=SPLACE
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir={#OutDir}
OutputBaseFilename=SPLACE-{#AppVersion}-windows-x64-setup
SetupIconFile={#Icon}
UninstallDisplayIcon={app}\splace.exe
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern

[Languages]
Name: "en"; MessagesFile: "compiler:Default.isl"
Name: "pt"; MessagesFile: "compiler:Languages\BrazilianPortuguese.isl"
Name: "es"; MessagesFile: "compiler:Languages\Spanish.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "{#Bundle}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\SPLACE"; Filename: "{app}\splace.exe"
Name: "{autodesktop}\SPLACE"; Filename: "{app}\splace.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\splace.exe"; Description: "{cm:LaunchProgram,SPLACE}"; Flags: nowait postinstall skipifsilent
