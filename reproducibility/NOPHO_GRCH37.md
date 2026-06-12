# NOPHO GRCh37/hg19 sensitivity reanalysis

The corrected script is `scripts/analyze_nopho_cnv_grch37.py`.

## Inputs

- 262 public CNVkit `.cns` files from Zenodo records `15167703`, `15173882`,
  and `15174016`.
- Standard public Boyle-Lab hg19 blacklist:
  `https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz`
- Blacklist SHA-256:
  `1a4ba636f791936ab8952cb068f496ccbd55ec4753539547c5d3d055ed00642a`
- Official study-code commit:
  `48607138ce0df2950d3ffc0f304246c38acf023f`

The exact study-specific `hg19-blacklist.v2.mod.bed` file was not distributed
with the public records. The analysis is therefore described as
official-workflow-compatible rather than an exact reproduction.

## Run

```powershell
python .\scripts\analyze_nopho_cnv_grch37.py `
  --source-root "PATH_TO_NOPHO_RAW_ROOT" `
  --blacklist "PATH_TO_hg19-blacklist.v2.bed.gz" `
  --output "results\nopho_grch37"
```

## Primary workflow-compatible definition

- deletion: `log2 <= -0.4`
- gain: `log2 >= 0.3`
- merge adjacent same-type events when gap is `<10 kb`
- retain events with size `>=1 Mb`
- discard events with at least 50% overlap with one blacklist interval
- use GRCh37 gene coordinates for overlap annotation

Gene overlap does not establish a breakpoint, minimal commonly deleted region,
or causal target.
