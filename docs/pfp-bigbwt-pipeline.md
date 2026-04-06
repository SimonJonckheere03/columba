# PFP and Big-BWT Pipeline

## Overview

The current Columba PFP build uses two different concerns:

- PFP remains the build mechanism that produces the phrase parse and the
  Big-BWT payload text.
- Runtime liftover is now defined directly on the plain searchable payload
  sequences, not on the internal sentinel-padded PFP text.

That means the runtime coordinate system is the concatenation of plain
per-sequence payload strings:

```text
forward_payload = seq0 || seq1 || seq2 || ...
```

Columba still requires both:

- the forward BWT built on `forward_payload`
- the reverse BWT built on `reverse_payload = reverse(forward_payload)`

## Runtime Coordinate Space

Runtime alignment now operates only in payload-space.

### Payload-space sequence inventory

The searchable payload is represented as a concatenation of named sequences.
Those sequence boundaries are described by:

- `<base>.pos`
- `<base>.sna`
- `<base>.fsid`
- `<base>.headerSN.bin`

Columba already uses these files in the plain multi-FASTA build, so the PFP
build now reuses the same mechanism.

### `.pos`

Binary array of sequence start offsets in the forward payload concatenation.

If the payload is:

```text
seq0 || seq1 || seq2
```

then `.pos` stores:

```text
[0, len(seq0), len(seq0)+len(seq1), len(seq0)+len(seq1)+len(seq2)]
```

At runtime, Columba uses `.pos` to decide:

- which payload sequence contains a hit
- whether a hit crosses a sequence boundary
- the sequence-relative offset of an accepted hit

### `.sna`

Length-prefixed sequence names in the same order as `.pos`.

### `.fsid`

The first sequence ID per logical source file/group.

For the current PFP-driven build this is written as one logical group, so the
runtime file contains a single value `0`.

### `.headerSN.bin`

The SAM `@SQ` header lines for the primary reporting coordinate system.

Because single-end liftover reports primary alignments in lifted reference
space, `.headerSN.bin` contains the lifted reference contig names and lifted
reference lengths, not the haplotype names.

## Liftover Files

### `.lidx`

`.lidx` is now a payload-space sequence manifest.

It stores:

- one entry per payload sequence
- the payload sequence name
- the payload sequence length

The order matches:

- `.pos`
- `.sna`
- `.ldx`
- the forward payload concatenation

### `.ldx`

`.ldx` is now payload-space LevioSAM metadata.

It stores:

- the payload-space sequence-start structure
- the payload-space sequence names
- the number of reference contigs
- one LevioSAM `lift::Lift` per payload sequence
- one reference-payload offset per payload sequence

For each payload sequence:

- reference sequences get an identity lift
- haplotype sequences get a LevioSAM lift built from the selected VCF alleles

Each lift maps:

```text
sequence-relative payload offset -> sequence-relative reference offset
```

The stored offset then converts that lifted coordinate into the global lifted
reference concatenation.

### `.bbwt.map`

`.bbwt.map` and `.rev.bbwt.map` may still be emitted as debugging sidecars by
the PFP build, but they are no longer a runtime dependency.

Columba does not load them in the payload-space redesign.

## Artifact Table

| Artifact | Producer | Coordinate space | Purpose | Current consumer |
| --- | --- | --- | --- | --- |
| `<base>.parse` / `.dict` / `.last` / `.sai` / `.occ` | `pfp++64` | internal PFP text | Forward PFP artifacts | Debug / parity |
| `<base>.rev.parse` / `.rev.dict` / `.rev.last` / `.rev.sai` / `.rev.occ` | `pfp++64` | internal reverse PFP text | Reverse internal/debug PFP artifacts | Debug |
| `<base>.lidx` | `pfp++64` | payload-space | Payload sequence manifest | Columba runtime |
| `<base>.ldx` | `pfp++64` | payload-space | Payload-space LevioSAM lifts | Columba runtime |
| `<base>.bbwt` | `pfp++64` | payload-space | Forward Big-BWT payload text | Big-BWT |
| `<base>.rev.bbwt` | `pfp++64` | reverse payload-space | Reverse Big-BWT payload text | Big-BWT |
| `<base>.pos` / `.sna` / `.fsid` | `columba_build_pfp.sh` | payload-space | Runtime sequence boundaries and names | Columba runtime |
| `<base>.headerSN.bin` | `columba_build_pfp.sh` | lifted reference space | SAM header lines | Columba runtime |
| `<base>.bbwt.parse` / `.dict` / `.last` / `.sai` / `.occ` | `newscanNT.x` | payload-space | Forward Big-BWT parser artifacts | `bwtparse64`, `pfbwtNT64.x` |
| `<base>.rev.bbwt.parse` / `.dict` / `.last` / `.sai` / `.occ` | `newscanNT.x` | reverse payload-space | Reverse Big-BWT parser artifacts | `bwtparse64`, `pfbwtNT64.x` |
| `<base>.bwt` / `.ssa` / `.esa` | Big-BWT handoff | payload-space | Final forward BWT and SA samples | `columba_build --pfp` |
| `<base>.rev.bwt` / `.rev.ssa` / `.rev.esa` | Big-BWT handoff | reverse payload-space | Final reverse BWT and SA samples | `columba_build --pfp` |
| `<base>.cct`, `.smpf`, `.smpl`, `.ftr`, `.ltr`, `.rev.smpf`, `.rev.smpl` | `columba_build --pfp` | derived from final BWTs | Columba runtime index data | Columba runtime |

## Build Pipeline

1. `pfp++64` reads the reference FASTA and VCF, materializes the payload
   sequences, builds the forward PFP artifacts, and writes payload-space
   `.lidx` / `.ldx`.
2. `pfp++64` also writes:
   - `<base>.bbwt`
   - `<base>.rev.bbwt`
3. `columba_build_pfp.sh` writes:
   - `.pos`
   - `.sna`
   - `.fsid`
   - `.headerSN.bin`
   directly from the payload-space sequence inventory and the reference FASTA.
4. `newscanNT.x`, `bwtparse64`, and `pfbwtNT64.x` build the forward BWT from
   `<base>.bbwt`.
5. The same tools build the reverse BWT from `<base>.rev.bbwt`.
6. `columba_build --pfp` builds the Columba runtime structures from the final
   forward and reverse BWTs.

## Runtime Single-End Liftover

For each accepted single-end hit:

1. search returns a hit in forward payload coordinates
2. Columba uses `.pos` to assign the hit to one payload sequence
3. if the hit crosses a payload-sequence boundary, Columba rejects it
4. Columba converts the hit to a sequence-relative payload offset
5. Columba loads the matching `lift::Lift` from `.ldx`
6. Columba lifts that payload offset directly into reference space
7. the primary SAM coordinate is reported in lifted reference space
8. the payload/haplotype origin is retained in the `XV:Z:` tag

If multiple haplotype hits lift to the same reference placement:

- one lifted primary alignment is kept
- all distinct origin placements are retained in `XV:Z:`

## Reverse Text Requirement

The reverse build is still mandatory.

Once the forward payload string is fixed, the reverse text is simply:

```text
reverse_payload = reverse(forward_payload)
```

Big-BWT then builds:

- `<base>.rev.bwt`
- `<base>.rev.ssa`
- `<base>.rev.esa`

The runtime redesign removes the old payload-to-PFP hop, but it does not
remove Columba’s need for the reverse searchable index.

## Move Structure Packing Note

The Columba move structures store an extra sentinel row and therefore need
bit widths sized for inclusive boundary values:

- `textSize + 1`
- `nrOfRuns + 1`

That fix remains required in:

- `MoveLFReprBP::initialize/load`
- `MovePhiReprBP::initialize/load`

Without it, exact power-of-two sizes truncate the sentinel boundary values and
corrupt the move representation.

## Current Scope

- Single-end liftover is supported.
- Paired-end liftover is intentionally unsupported in this version.
- Runtime does not depend on `.bbwt.map`.
