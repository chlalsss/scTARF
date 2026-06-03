$ErrorActionPreference = "Continue"

$Root = "D:\RStudio\work\module4_traitaware_pseudopair"
$Trait = "MCH"
$Tag = "v3_sparse_cite_erythroid_raw_atac_full_pipeline_dropout"
$RunLabel = "cite_erythroid"
$BatchSize = "64"
$LogDir = Join-Path $Root "05_robustness\logs"

New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

Set-Location $Root

function Test-ConditionComplete {
    param([double]$PDrop, [int]$Seed)

    $runId = "p{0:0.00}_seed{1}" -f $PDrop, $Seed
    $csv = Join-Path $Root ("05_robustness\output\{0}\runs\{1}\04_evaluation\{2}_neighborhood_fdr_ge0p80\{1}_v3_sparse_{2}_full_pipeline_auroc.csv" -f $Tag, $runId, $RunLabel)
    if (-not (Test-Path -LiteralPath $csv)) {
        return $false
    }

    $hits = Select-String -LiteralPath $csv -Pattern '(^|,)"?MCH"?(,|$)' -Quiet
    return [bool]$hits
}

function Invoke-OneCondition {
    param([double]$PDrop, [int]$Seed)

    $runId = "p{0:0.00}_seed{1}" -f $PDrop, $Seed
    $log = Join-Path $LogDir ("module4_MCH_{0}_{1}.log" -f $runId, (Get-Date -Format "yyyyMMdd_HHmmss"))

    "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] START $runId" | Tee-Object -FilePath $log
    $env:MODULE4_TRAIT = $Trait
    $env:MODULE4_P_DROP = [string]$PDrop
    $env:MODULE4_SEED = [string]$Seed
    $env:MODULE4_BATCH_SIZE = $BatchSize

    & Rscript _run_module4_one_condition.R 2>&1 | Tee-Object -FilePath $log -Append
    $code = $LASTEXITCODE

    Remove-Item Env:\MODULE4_TRAIT -ErrorAction SilentlyContinue
    Remove-Item Env:\MODULE4_P_DROP -ErrorAction SilentlyContinue
    Remove-Item Env:\MODULE4_SEED -ErrorAction SilentlyContinue
    Remove-Item Env:\MODULE4_BATCH_SIZE -ErrorAction SilentlyContinue

    if ($code -ne 0) {
        "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] FAIL $runId exit=$code" | Tee-Object -FilePath $log -Append
        exit $code
    }

    if (-not (Test-ConditionComplete -PDrop $PDrop -Seed $Seed)) {
        "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] FAIL $runId missing expected AUROC CSV row" | Tee-Object -FilePath $log -Append
        exit 2
    }

    "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] DONE $runId" | Tee-Object -FilePath $log -Append
}

"[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] MCH missing-condition runner started"

$pDrops = 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
$seeds = 1,2,3

foreach ($p in $pDrops) {
    foreach ($seed in $seeds) {
        $runId = "p{0:0.00}_seed{1}" -f $p, $seed
        if (Test-ConditionComplete -PDrop $p -Seed $seed) {
            "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] SKIP $runId already complete"
            continue
        }
        Invoke-OneCondition -PDrop $p -Seed $seed
    }
}

"[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] MCH missing-condition runner finished"
