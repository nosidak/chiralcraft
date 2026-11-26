param(
  [Parameter(Mandatory=$true)][string]$Version,
  [Parameter(Mandatory=$true)][ValidateSet("x64","arm64")][string]$Arch,
  [Parameter(Mandatory=$true)][string]$BinPath,
  [string]$OutDir = "dist"
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$BinFull = (Resolve-Path $BinPath).Path
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$OutFull = (Resolve-Path $OutDir).Path

$PkgBase = "chiralcraft-$Version-windows-$Arch"

# 1. Portable ZIP
Write-Host "==> Building $PkgBase-portable.zip"
$ZipStage = Join-Path $env:TEMP "$PkgBase-portable"
if (Test-Path $ZipStage) { Remove-Item -Recurse -Force $ZipStage }
New-Item -ItemType Directory -Force -Path $ZipStage | Out-Null
Copy-Item $BinFull (Join-Path $ZipStage "chiralcraft.exe")
Copy-Item (Join-Path $RepoRoot "README.md") $ZipStage
Copy-Item (Join-Path $RepoRoot "LICENSE") $ZipStage
Copy-Item (Join-Path $RepoRoot "CHANGELOG.md") $ZipStage
$ZipOut = Join-Path $OutFull "$PkgBase-portable.zip"
if (Test-Path $ZipOut) { Remove-Item -Force $ZipOut }
Compress-Archive -Path "$ZipStage\*" -DestinationPath $ZipOut

# 2. WiX MSI Installer
$WixArch = if ($Arch -eq "arm64") { "arm64" } else { "x64" }
$WixCmd = Get-Command wix -ErrorAction SilentlyContinue
if ($WixCmd) {
  Write-Host "==> Building $PkgBase.msi"
  $MsiOut = Join-Path $OutFull "$PkgBase.msi"
  & wix build (Join-Path $PSScriptRoot "chiralcraft.wxs") `
    -arch $WixArch `
    -d "ProductVersion=$Version" `
    -d "BinPath=$BinFull" `
    -d "RepoRoot=$RepoRoot" `
    -o $MsiOut
  if ($LASTEXITCODE -ne 0) { throw "wix build failed with exit code $LASTEXITCODE" }
}

# 3. Inno Setup EXE Installer
$Iscc = Get-Command iscc -ErrorAction SilentlyContinue
if (-not $Iscc) {
  $Candidate = "C:\Program Files (x86)\Inno Setup 6\ISCC.exe"
  if (Test-Path $Candidate) { $Iscc = $Candidate }
}
if ($Iscc) {
  Write-Host "==> Building $PkgBase-setup.exe"
  $InnoAllowed = if ($Arch -eq "arm64") { "arm64" } else { "x64compatible" }
  & $Iscc `
    "/DMyAppVersion=$Version" `
    "/DMyAppArch=$Arch" `
    "/DInnoArchAllowed=$InnoAllowed" `
    "/DBinPath=$BinFull" `
    "/DRepoRoot=$RepoRoot" `
    "/DOutDir=$OutFull" `
    (Join-Path $PSScriptRoot "chiralcraft.iss")
  if ($LASTEXITCODE -ne 0) { throw "ISCC failed with exit code $LASTEXITCODE" }
}

Get-ChildItem $OutFull -Filter "$PkgBase*" | Format-Table Name, Length
