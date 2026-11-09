# NJoiner

NJoiner reconstructs non-homologous end-joining (NHEJ) repair junctions from their DNA sequences. Given the original sequence with the double-strand break marked, and the sequence of a repair product, it aligns the product to each side of the break and reports the microhomology, inserted (filler) bases and deleted bases at the junction, and whether any filler has a template nearby. It was written for classifying large collections of junctions objectively and repeatably, so that collections from different experimental conditions can be compared.

NJoiner was written by William R. Engels at the University of Wisconsin–Madison in 2010 and presented with Carlos C. Flores at a Keystone Symposium in 2011 (poster included). It is the source of the junction classes in Table 3 of *Genomic Position Effects on Double-Strand Break Repair in Drosophila* (Preston, Flores, Johnson-Schlitz, Hum and Engels; unpublished manuscript, 2013/2026).

The original program is a Mac OS X 10.6 Cocoa application that no longer runs on current macOS. This repository preserves its source and provides a Python port of the alignment and classification logic that reproduces the original program's output.

## Quick start

Requires only Python 3 (standard library).

```
python3 python/njoiner.py test/No109Rr3EJ.txt -o out.csv --html out.html
```

Input format: a line `Original, <sequence>` with the cut marked `^` (blunt) or `^XXXX^` (staggered, overhang between the carets), then one line per repair product, `name, sequence`. Parameter sets: `--params poster2011` (default), `april2010`, `output2010`.

To check the port against the original program's 2010 output and re-derive Table 3 of the manuscript:

```
python3 python/validate.py
```

## Contents

| Path | Contents |
|:--|:--|
| `NJoiner - description.md` | The method, its parameter settings, known quirks of the original, and how the port was checked |
| `2011 Engels & Flores - NJoiner poster (Keystone).pdf` | The authors' poster, Keystone Symposium, January–February 2011 |
| `source (Objective-C, 2010)/` | Original Cocoa source and Xcode project. The alignment logic is in `NJAligner.m` and `NJoin.m`. |
| `python/njoiner.py` | Python 3 port of the alignment core, with text and colour HTML displays |
| `python/validate.py` | Checks the port against the 2010 output and re-derives Table 3 |
| `test/` | Test input, NJoiner's own output for it (12 April 2010), the port's output, and the comparison |
| `table3/` | The port's re-derivation of Table 3: per-junction classes, alignments and comparison with the authors' tabulation |
| `data/NHEJ_junction_sequences/` | Junction sequences from the *Rr3* reporter in *Drosophila*, as NJoiner input files |

## Algorithm references

Smith, T. F., and M. S. Waterman, 1981 Identification of common molecular subsequences. *J. Mol. Biol.* 147: 195–197.

Gotoh, O., 1982 An improved algorithm for matching biological sequences. *J. Mol. Biol.* 162: 705–708.

## Licence

MIT; see `LICENSE`.

This repository was assembled in 2026 by William R. Engels with the assistance of Claude, an AI model made by Anthropic, which wrote the Python port and the documentation; W. R. Engels reviewed them.
