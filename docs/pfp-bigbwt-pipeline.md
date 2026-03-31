# PFP and Big-BWT Pipeline

## Overview

The current Columba PFP build has two separate text representations:

- `pfp_pos`: the full forward or reverse PFP text, including structural sentinel bytes.
- `bbwt_pos`: the stripped Big-BWT payload text, containing only biological characters.

The forward PFP coordinate space is the source of truth for:

- forward `.parse`, `.dict`, `.last`, `.sai`, `.occ`
- forward `.lidx`
- forward `.ldx`

The Big-BWT-facing payloads are derived views:

- `<base>.bbwt`
- `<base>.rev.bbwt`

These payloads are built by removing structural PFP sentinels from the fully materialized forward and reverse texts. Big-BWT indexes the payloads, not the sentinel-rich PFP texts.

## Coordinate Spaces

### `pfp_pos`

`pfp_pos` is the coordinate system of the exact materialized PFP text. It includes:

- `DOLLAR` (`0x02`)
- `DOLLAR_SEQUENCE` (`0x04`)
- `DOLLAR_PRIME` (`0x05`)

This is the space used by the forward `.lidx` and `.ldx` outputs.

### `bbwt_pos`

`bbwt_pos` is the coordinate system of the stripped payload passed to Big-BWT. It excludes all structural sentinels and contains only biological characters after current preprocessing rules.

### Mapping

Each payload has a matching packed keep-bitvector:

- `<base>.bbwt.map`
- `<base>.rev.bbwt.map`

File format:

1. `uint64_t full_pfp_length`
2. `uint64_t kept_count`
3. packed bitvector bytes over `full_pfp_length`

Bit meaning:

- `1`: the character at that `pfp_pos` was kept in the Big-BWT payload
- `0`: the character was dropped because it is a structural sentinel

Expected load-time operations:

- `bbwt_pos -> pfp_pos`: `select1(bbwt_pos + 1)`
- `pfp_pos -> bbwt_pos`: `rank1(pfp_pos + 1) - 1` if the bit is `1`

Boundary validation rule:

- map the start and end of an interval back to `pfp_pos`
- reject the interval if the mapped span is longer than the payload interval
- this detects removed sentinel bytes inside the candidate interval

## Artifact Table

| Artifact | Producer | Coordinate space | Purpose | Current consumer |
| --- | --- | --- | --- | --- |
| `<base>.parse` | `pfp++64` | `pfp_pos` | Forward phrase sequence | Debug/parity only |
| `<base>.dict` | `pfp++64` | `pfp_pos` phrases | Forward phrase dictionary | Debug/parity only |
| `<base>.last` / `<base>.sai` / `<base>.occ` | `pfp++64` | `pfp_pos` parse stream | Forward PFP metadata | Debug/parity only |
| `<base>.rev.parse` / `.rev.dict` / `.rev.last` / `.rev.sai` / `.rev.occ` | `pfp++64` | reverse `pfp_pos` | Reverse internal/debug PFP artifacts | Debug only |
| `<base>.lidx` | `pfp++64` | forward `pfp_pos` | Sequence index lengths/names | Future alignment/liftover |
| `<base>.ldx` | `pfp++64` | forward `pfp_pos` | Lift index | Future alignment/liftover |
| `<base>.bbwt` | `pfp++64` | `bbwt_pos` | Forward stripped payload for Big-BWT | Big-BWT |
| `<base>.bbwt.map` | `pfp++64` | `pfp_pos` bitvector | Forward `bbwt_pos -> pfp_pos` transform | Future alignment/liftover |
| `<base>.rev.bbwt` | `pfp++64` | reverse `bbwt_pos` | Reverse stripped payload for Big-BWT | Big-BWT |
| `<base>.rev.bbwt.map` | `pfp++64` | reverse `pfp_pos` bitvector | Reverse `bbwt_pos -> pfp_pos` transform | Future alignment/liftover |
| `<base>.pos` / `.sna` / `.fsid` / `.headerSN.bin` | `columba_build_pfp.sh` | `bbwt_pos` | Runtime sequence names, starts, file grouping, and SAM header lines derived from `.lidx` + `.bbwt.map` | Columba runtime |
| `<base>.bbwt.parse` / `.dict` / `.last` / `.sai` / `.occ` | `newscanNT.x` | `bbwt_pos` | Forward Big-BWT parser artifacts | `bwtparse64`, `pfbwtNT64.x` |
| `<base>.rev.bbwt.parse` / `.dict` / `.last` / `.sai` / `.occ` | `newscanNT.x` | reverse `bbwt_pos` | Reverse Big-BWT parser artifacts | `bwtparse64`, `pfbwtNT64.x` |
| `<base>.bwt` / `.ssa` / `.esa` | Big-BWT handoff | `bbwt_pos` | Final forward BWT and SA samples | `columba_build --pfp` |
| `<base>.rev.bwt` / `.rev.ssa` / `.rev.esa` | Big-BWT handoff | reverse `bbwt_pos` | Final reverse BWT and SA samples | `columba_build --pfp` |
| `<base>.cct`, `.smpf`, `.smpl`, `.ftr`, `.ltr`, `.rev.smpf`, `.rev.smpl` | `columba_build --pfp` | derived from final BWTs | Columba runtime index data | Columba runtime |

## Liftover Semantics

The forward liftover structures remain sentinel-aware and unchanged.

### `.lidx`

`.lidx` stores contig names and contig lengths in the forward PFP coordinate system. The reported lengths include the PFP tail convention used during forward parsing.

### `.ldx`

`.ldx` stores the merged LevioSAM lift structures for the forward PFP coordinate system. It is built from the current forward parser path and must not be rewritten to match stripped Big-BWT payload coordinates.

### Consequence

Any future alignment code that obtains text positions from the final `.bwt` must:

1. interpret them first as `bbwt_pos`
2. map them back to `pfp_pos` with `.bbwt.map`
3. only then apply `.lidx` / `.ldx`

## Current Build Pipeline

1. `pfp++64` builds forward PFP artifacts, reverse debug artifacts, `.lidx`, `.ldx`, and the stripped `.bbwt` payloads plus `.bbwt.map`.
2. `columba_build_pfp.sh` derives `.pos`, `.sna`, `.fsid`, and `.headerSN.bin` from `.lidx` plus `.bbwt.map`, so the aligner sees the same sequence metadata contract as in the plain FASTA build.
3. `newscanNT.x` runs on `<base>.bbwt` and `<base>.rev.bbwt`.
4. `bwtparse64` and `pfbwtNT64.x` build payload-based Big-BWT outputs.
5. The final payload-based `.bwt/.ssa/.esa` files are installed to the standard names expected by `columba_build --pfp`.
6. `columba_build --pfp` builds the Columba runtime structures from the final BWT files.

## Move Structure Packing Note

The Columba move structures store one extra sentinel row in both the LF move
representation and the Phi move representation.

That sentinel row stores inclusive boundary values:

- `inputStartPos = textSize`
- `outputStartPos = textSize`
- `outputStartRun = nrOfRuns`

This means the bit widths in `src/bmove/moverepr.cpp` must be sized for
`textSize + 1` and `nrOfRuns + 1`, not just for `textSize` and `nrOfRuns`.

The current Columba version carries that fix in:

- `MoveLFReprBP::initialize/load`
- `MovePhiReprBP::initialize/load`

Without that fix, power-of-two cases are truncated during packing. For example,
`textSize = 32` would use only 5 bits and the stored sentinel boundary `32`
would wrap to `0`. In the PFP build path this caused the last run boundary to
be lost and later triggered a PLCP select assertion during
`constructRunLengthEncodedPLCP(...)`.

This fix is independent of PFP itself, but it is required for the PFP-based
index build to remain consistent with the Big-BWT-generated `.bwt/.ssa/.esa`
files.

## Downstream Expectations

Current Columba code consumes:

- `.bwt`
- `.rev.bwt`
- `.ssa`
- `.esa`
- `.rev.ssa`
- `.rev.esa`
- `.pos`
- `.sna`
- `.fsid`
- `.headerSN.bin`
- the derived Columba build outputs such as `.cct`, `.smpf`, `.smpl`, `.ftr`, `.ltr`

Current Columba code does **not** yet load `.lidx` or `.ldx`.

Expected future alignment flow:

1. locate a candidate interval in the forward or reverse BWT
2. convert the resulting text positions from `bbwt_pos` to `pfp_pos` with `.bbwt.map`
3. reject intervals whose mapped span crosses removed sentinels
4. apply sequence-boundary validation in `pfp_pos`
5. map accepted forward positions through `.ldx`

## Invariants

- Forward `.parse/.dict/.last/.sai/.occ/.lidx/.ldx` remain byte-identical to `moni-align`.
- `.bbwt` and `.rev.bbwt` contain no structural PFP sentinel bytes.
- `.bbwt.map` and `.rev.bbwt.map` reconstruct the kept payload positions exactly.
- Any interval that crosses a removed structural sentinel must be rejected after mapping back to `pfp_pos`.
- Move-structure sentinel rows must preserve `textSize` and `nrOfRuns` exactly; they may not be truncated during bit packing.
