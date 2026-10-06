param(
    [string]$Destination = (Join-Path $env:LOCALAPPDATA 'Programs/fast-gfn1-xtb'),
    [switch]$AddToPath
)
$ErrorActionPreference = 'Stop'
$source = [IO.Path]::GetFullPath($PSScriptRoot)
$target = [IO.Path]::GetFullPath($Destination).TrimEnd('\', '/')
if ($target.Equals($source, [StringComparison]::OrdinalIgnoreCase) -or
    $target.StartsWith($source + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) {
    throw 'Choose an installation folder outside the extracted package.'
}
if (Test-Path -LiteralPath $target) {
    throw "Installation folder already exists: $target. Choose another destination."
}
$manifest = Get-Content -LiteralPath (Join-Path $source 'PACKAGE_MANIFEST.json') -Raw -Encoding UTF8 | ConvertFrom-Json
$files = @()
foreach ($entry in $manifest.files.PSObject.Properties) {
    $path = [IO.Path]::GetFullPath((Join-Path $source $entry.Name))
    if (!$path.StartsWith($source + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) {
        throw 'Invalid package file path.'
    }
    if ((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash -ne $entry.Value) {
        throw "Package checksum failed: $($entry.Name)"
    }
    $files += $entry.Name
}
New-Item -ItemType Directory -Path $target | Out-Null
foreach ($name in @($files) + @('PACKAGE_MANIFEST.json')) {
    $destinationFile = Join-Path $target $name
    New-Item -ItemType Directory -Force -Path ([IO.Path]::GetDirectoryName($destinationFile)) | Out-Null
    Copy-Item -LiteralPath (Join-Path $source $name) -Destination $destinationFile
}
# Verify the installed files before changing the user's command search path.
foreach ($entry in $manifest.files.PSObject.Properties) {
    if ((Get-FileHash -LiteralPath (Join-Path $target $entry.Name) -Algorithm SHA256).Hash -ne $entry.Value) {
        throw "Installed checksum failed: $($entry.Name)"
    }
}
if ($AddToPath) {
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    $parts = @($userPath -split ';' | Where-Object { $_ })
    if (!($parts | Where-Object { $_.TrimEnd('\', '/') -ieq $target })) {
        [Environment]::SetEnvironmentVariable('Path', (($parts + @($target)) -join ';'), 'User')
    }
}
Write-Output "Installed: $target"
Write-Output 'Open a new terminal, then run: fast-gfn1-xtb input.xyz --gfn 1 --grad'
