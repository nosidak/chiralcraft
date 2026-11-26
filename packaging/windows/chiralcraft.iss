#ifndef MyAppVersion
  #define MyAppVersion "1.0.0"
#endif
#ifndef MyAppArch
  #define MyAppArch "x64"
#endif

[Setup]
AppId={{94E7A5D1-7C3B-4A68-9E8F-3C1A4D2B7E55}
AppName=ChiralCraft
AppVersion={#MyAppVersion}
AppPublisher=nosidak
AppPublisherURL=https://github.com/nosidak/chiralcraft
AppSupportURL=https://github.com/nosidak/chiralcraft/issues
DefaultDirName={autopf}\ChiralCraft
DefaultGroupName=ChiralCraft
LicenseFile={#RepoRoot}\LICENSE
OutputDir={#OutDir}
OutputBaseFilename=chiralcraft-{#MyAppVersion}-windows-{#MyAppArch}-setup
Compression=lzma2/ultra64
SolidCompression=yes
ArchitecturesAllowed={#InnoArchAllowed}
ArchitecturesInstallIn64BitMode={#InnoArchAllowed}
PrivilegesRequiredOverridesAllowed=dialog
WizardStyle=modern

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "{#BinPath}"; DestDir: "{app}"; DestName: "chiralcraft.exe"; Flags: ignoreversion
Source: "{#RepoRoot}\README.md"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#RepoRoot}\LICENSE"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#RepoRoot}\CHANGELOG.md"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\ChiralCraft"; Filename: "{app}\chiralcraft.exe"
Name: "{group}\Uninstall ChiralCraft"; Filename: "{uninstallexe}"
Name: "{autodesktop}\ChiralCraft"; Filename: "{app}\chiralcraft.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\chiralcraft.exe"; Description: "Launch ChiralCraft"; Flags: nowait postinstall skipifsilent
