param(
    [Parameter(Mandatory = $true)]
    [string]$SourceRoot,
    [string]$WorkspaceRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path,
    [int]$Permutations = 200,
    [string]$NophoRawRoot,
    [string]$Hg19Blacklist,
    [switch]$BuildManuscript
)

$ErrorActionPreference = "Stop"

python (Join-Path $WorkspaceRoot "scripts\run_scientific_reports_robustness_ml.py") `
    --workspace-root $WorkspaceRoot `
    --source-root $SourceRoot `
    --permutations $Permutations

if ($NophoRawRoot -and $Hg19Blacklist) {
    python (Join-Path $WorkspaceRoot "scripts\analyze_nopho_cnv_grch37.py") `
        --source-root $NophoRawRoot `
        --blacklist $Hg19Blacklist `
        --output (Join-Path $WorkspaceRoot "results\nopho_grch37")
}

if ($BuildManuscript) {
    python (Join-Path $WorkspaceRoot "scripts\generate_historical_only_single_center_figure.py") `
        --workspace-root $WorkspaceRoot

    python (Join-Path $WorkspaceRoot "scripts\integrate_scientific_reports_hard_fixes_v7.py") `
        --workspace-root $WorkspaceRoot
}
