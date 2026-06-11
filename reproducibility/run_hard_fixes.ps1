param(
    [Parameter(Mandatory = $true)]
    [string]$SourceRoot,
    [string]$WorkspaceRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [int]$Permutations = 200,
    [switch]$BuildManuscript
)

$ErrorActionPreference = "Stop"

python (Join-Path $WorkspaceRoot "scripts\run_scientific_reports_robustness_ml.py") `
    --workspace-root $WorkspaceRoot `
    --source-root $SourceRoot `
    --permutations $Permutations

if ($BuildManuscript) {
    python (Join-Path $WorkspaceRoot "scripts\generate_historical_only_single_center_figure.py") `
        --workspace-root $WorkspaceRoot

    python (Join-Path $WorkspaceRoot "scripts\integrate_scientific_reports_hard_fixes_v7.py") `
        --workspace-root $WorkspaceRoot
}
