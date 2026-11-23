; OpenProblemsLab for Windows - per-user installer (no admin rights). build.py compiles it with:
;   iscc /DAppVersion=1.0.0 /DSourceDir=<dist\OpenProblemsLab> /DOutputDir=<dist> /DIconFile=<build\openproblemslab.ico> openproblemslab.iss
; Installs to %LOCALAPPDATA%\Programs\OpenProblemsLab, adds OpenProblemsLab to the Start menu (and optionally the desktop) and lists it
; in Settings > Apps. The window keeps its data (what you saved in it) in %LOCALAPPDATA%\OpenProblemsLab, which
; uninstalling leaves in place so a reinstall finds it again.

#ifndef AppVersion
  #define AppVersion "1.0.0"
#endif
#ifndef SourceDir
  #define SourceDir "dist\OpenProblemsLab"
#endif
#ifndef OutputDir
  #define OutputDir "dist"
#endif
#ifndef IconFile
  #define IconFile "build\openproblemslab.ico"
#endif

[Setup]
AppId={{E81247BA-698B-569E-A6DE-E1E1C366EEFB}
AppName=OpenProblemsLab
AppVersion={#AppVersion}
AppVerName=OpenProblemsLab {#AppVersion}
AppPublisher=Juwen Lie
AppPublisherURL=https://juw12.github.io/open-problems-lab/
AppSupportURL=https://github.com/juw12/open-problems-lab
DefaultDirName={autopf}\OpenProblemsLab
DisableDirPage=yes
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir={#OutputDir}
OutputBaseFilename=OpenProblemsLab-Windows-Setup
SetupIconFile={#IconFile}
UninstallDisplayIcon={app}\OpenProblemsLab.exe
UninstallDisplayName=OpenProblemsLab
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
CloseApplications=no
SetupLogging=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{userprograms}\OpenProblemsLab"; Filename: "{app}\OpenProblemsLab.exe"; Comment: "Reproducible computations on open problems in number theory"
Name: "{userdesktop}\OpenProblemsLab"; Filename: "{app}\OpenProblemsLab.exe"; Comment: "Reproducible computations on open problems in number theory"; Tasks: desktopicon

[Run]
Filename: "{app}\OpenProblemsLab.exe"; Description: "{cm:LaunchProgram,OpenProblemsLab}"; Flags: nowait postinstall skipifsilent

[UninstallRun]
Filename: "{app}\OpenProblemsLab.exe"; Parameters: "--stop"; Flags: runhidden waituntilterminated; RunOnceId: "StopOpenProblemsLab"

[Code]
function PrepareToInstall(var NeedsRestart: Boolean): String;
var
  ResultCode: Integer;
begin
  Result := '';
  { Updating: close an open OpenProblemsLab window first so its files can be replaced. }
  if FileExists(ExpandConstant('{app}\OpenProblemsLab.exe')) then
    Exec(ExpandConstant('{app}\OpenProblemsLab.exe'), '--stop', '', SW_HIDE, ewWaitUntilTerminated, ResultCode);
end;
