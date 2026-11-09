---
title: "NJoiner: a program for reconstructing and classifying end-joining junctions"
author: William R. Engels (program, 2010); description written 2026
done: true
---

# NJoiner — description

*About this note (2026).* This description was written in 2026 by William R. Engels with the assistance of Claude, an AI model made by Anthropic. Claude read the 2010 Objective-C source in `source (Objective-C, 2010)/`, the authors' 2011 poster, and the program's surviving output; wrote the Python port in `python/`; and drafted this note. W. R. Engels reviewed it and takes responsibility for it.

## What NJoiner is

NJoiner was written by W. R. Engels in 2010 to classify large numbers of non-homologous end-joining (NHEJ) repair products objectively and repeatably. It was used in his laboratory at the University of Wisconsin–Madison on junctions from the *Rr3* reporter, in which I-*Sce*I makes a break in the *Drosophila* male germline and the repair products are recovered in the next generation and sequenced. It was presented as a poster by William Engels and Carlos Flores at a Keystone Symposium (30 January – 4 February 2011), included in this folder. The poster reports that the settings were tuned on a collection of 1337 end-joined sequences from the laboratory and from published work. In the Position Effects manuscript it is the source of the junction classes summarized in Table 3.

The motivating problem, as the poster puts it, is that the same repair product can often be aligned to the original sequence in more than one plausible way, and the choice changes whether it is scored as having microhomology, how long the microhomology is, and how many bases were lost on each side. No single event can be assigned with certainty; what matters is that collections of events from different conditions are classified by the same rule, so that they can be compared.

The program is a Cocoa application written in Objective-C for Mac OS X 10.6 (the poster says 10.5 or later). It was not written in Swift, which did not appear until 2014.

## Input

1. The original (substrate) sequence, with the cut marked. A blunt cut is written `^`; a staggered cut is written `^XXXX^` with the overhang between the carets. For *Rr3* the original is
   `…CTAGGG^ATAA^CAGGG…`, the 4-base 3′ overhang ATAA left by I-*Sce*I.
2. One or more repair-product sequences, from a text file (one `name, sequence` line per product), the clipboard, or a folder of ABI `.ab1` chromatogram files.

## The algorithm

Following the poster and the source (`NJAligner.m`, 10 September 2010; `NJoin.m`, 7 September 2010):

1. The original is split at the break into a left part and a right part. Any overhang is included in both parts.
2. The product is aligned separately against the left part and against the right part, using Smith–Waterman–Gotoh local alignment (Smith and Waterman 1981; Gotoh 1982) with affine gap penalties.
3. Product bases covered by both local alignments are called **microhomology**. Bases of the left part lying to the right of its alignment, and bases of the right part lying to the left of its alignment, are **missing** (deleted). Product bases lying between the two alignments are **filler** (inserted bases).
4. The two alignments are merged into a single three-row reconstruction (left end / product / right end), introducing gaps where needed so that every matched base of either alignment is kept.
5. If the filler is at least a minimum length (6 bases by default), it is aligned against the original sequence in both orientations to look for a template. A match scoring at least 70% of the maximum possible is reported as a direct (1) or inverted (2) template; otherwise 0.
6. The reconstruction is written out as a text diagram and as a row of a table.

## Output fields

| Field | Meaning |
|:--|:--|
| complex | 1 if the reconstruction needs internal gaps, or if matching resumes on one side after a stretch of non-matching bases |
| microhomology, micro length | the bases shared by both alignments, and their number |
| imperfect micro | intended to flag microhomology containing a mismatch or gap (but see Quirks) |
| filler, filler length | inserted bases between the two alignments |
| filler template | 0 none, 1 direct, 2 inverted (step 5) |
| missing left, missing right | bases deleted from each side of the break |

In the Position Effects manuscript, a junction with filler is classed as an **in-del** and one without filler as an **EJ** (a pure deletion). In NJoiner's output no junction has both filler and microhomology, so microhomology is scored only within the EJ class.

## Parameter settings

Three sets survive, and they differ:

| | Match | Mismatch | Ambiguous (N) | Deletion open | Deletion extend | Insertion open | Insertion extend |
|:--|--:|--:|--:|--:|--:|--:|--:|
| Saved in `parameters.plist` and `userPrefs.flist` (12–14 April 2010) | 4 | −8 | 0 | −16 | −2 | −12 | −2 |
| Reproduces the program's own output of 12 April 2010 | 4 | −8 | 0 | **−14** | −2 | −12 | −2 |
| Poster default (January 2011) | 4 | −8 | 0 | −16 | −4 | −12 | −4 |

All three use a minimum filler length of 6 for the template search and a 70% match threshold. The second row was found by trial: with the saved settings the port disagrees with the 2010 output on 12 of its 69 junctions (as it does with the poster settings), and with deletion-open −14 it agrees on all of them. The settings in force when the output was written were evidently not those later saved. The poster set is the one closest to the tabulation behind Table 3, and it is the port's default (`--params poster2011`).

## Quirks of the 2010 program

These were kept deliberately in the Python port, so that it reproduces what NJoiner did rather than what it was meant to do.

- **The ambiguity test uses the wrong index.** In the alignment loop the test for an N reads `a[i]=='N' || b[i]=='N'`, where `b` should be indexed by `j`. The reference never contains N, so in practice only an N in the product is treated as ambiguous; past the end of `b` the original read unallocated memory.
- **Uninitialized backtrace pointer.** The variable `pointer` that steers the backtrace is not initialized. The starting cell of a local maximum always has a diagonal backtrace, so the port starts it at 0.
- **Filler location in mixed coordinates.** The filler start is taken as the smaller of a left-part index and a product index.
- **The "imperfect micro" column is inverted.** In the 2010 output it is 1 when the junction has microhomology and that microhomology is *perfect*, and 0 otherwise. Its values should not be read as the label says.

## Status of the original program

The source is complete, including the interface files (`.nib`, `.xib`) and the Xcode project. The project is set to build against the Mac OS X 10.6 SDK with Objective-C garbage collection *required*. Garbage-collected applications are no longer supported by macOS, so the program will not run on a current Mac without rewriting its memory management, and its interface would also need rebuilding. Neither is necessary to use the method: the alignment and classification logic is in two files (`NJAligner.m`, `NJoin.m`) and has been ported.

## The Python port

`python/njoiner.py` is a statement-by-statement port of the alignment core. It does not reproduce the interface, the `.ab1` reader or the multithreading.

```
python3 python/njoiner.py INPUT.txt -o OUT.csv --html OUT.html --diagrams OUT.txt [--params poster2011|april2010|output2010]
```

`INPUT.txt` is an NJoiner input file: a line `Original, <sequence with ^ marks>` (or `LEFT`/`RIGHT`/`OVERHANG` lines), then one line per product, `name, sequence`. Lines whose name begins with `//` are comments. The text diagram is the five-line display of the original program (left end / match line / product / match line / right end); the HTML version colours missing bases, microhomology, filler and the overhang.

`python/validate.py` (run from this folder) checks the port in two ways.

1. **Against the program's own output.** For the test file `test/No109Rr3EJ.txt`, the port with deletion-open −14 matches NJoiner's output of 12 April 2010 in every field for all 69 junctions, except the "imperfect micro" column, which matches a perfect-microhomology flag in all 69 (see Quirks). The comparison is in `test/No109Rr3EJ comparison.csv`.
2. **Against Table 3.** Using the plate files in `data/NHEJ_junction_sequences/` and the poster parameters, it re-derives the classes behind Table 3 and compares them with the authors' tabulation. The results are in `table3/` (per-junction classes and HTML alignments for each sample).

| Sample | Junctions, tabulation / port | % of EJ with mh, tabulation / port |
|:--|--:|--:|
| 21C, young | 38 / 37 | 82.8 / 85.7 |
| 21C, old | 20 / 18 | 50.0 / 46.2 |
| 82D, young | 13 / 12 | 33.3 / 28.6 |
| 48C, young | 91 / 91 | 82.5 / 82.8 |
| 48C, 3 weeks | 89 / 87 | 76.7 / 76.7 |
| 48C, 4–6 weeks | 82 / 80 | 61.0 / 65.5 |
| 48C, held to old age | 78 / 67 | 45.5 / 65.0 |

The agreement is close except for the oldest 48C sample, whose surviving sequence file holds 11 fewer records than the tabulation; the missing records cannot be identified, so that row of the port should not be used. The site differences and the age trend at 48C up to 4–6 weeks are reproduced.

## Findings from the re-analysis

- **The 4-base microhomology class is not the overhang.** Although four bases is the length of the I-*Sce*I overhang, 28 of the 33 four-base joins at 21C and 48C (young males) pair the AGGG repeat present on both sides of the cut (`…CT`**`AGGG`**`^ATAA^C`**`AGGG`**`TAAT…`). These junctions lack the nine bases ATAACAGGG, so the ends must be trimmed back past the overhang before they can pair. Only the single four-base join at 82D used the ATAA overhang itself.
- **AGGG joins at 82D.** None of the 7 end joins at 82D used the AGGG repeat, against 28 of 90 at 21C and 48C combined (Fisher's exact test, two-tailed *P* = 0.19). The sample is too small to say which step of end joining is affected at 82D.
- **A probable duplicate.** Three of the 82D in-del reads (82D.a.9a, 82D.a.24b, 82D.a.44a) are identical, including the positions of their N bases, which suggests one sequence was entered three times.

## Files

See `README.md` in this folder.

## References

Gotoh, O., 1982 An improved algorithm for matching biological sequences. *J. Mol. Biol.* 162: 705–708.

Smith, T. F., and M. S. Waterman, 1981 Identification of common molecular subsequences. *J. Mol. Biol.* 147: 195–197.
