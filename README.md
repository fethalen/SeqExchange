SeqExchange
-----------

## About

SeqExchange is a Python script for assigning sequences from one species to
another. It is meant as a tool for simulating cross-contamination and an
abnormal number of paralogs (by paralogs, we here mean 2+ sequences belonging
to 1 species in an alignment).

SeqExchange has two modes: one that simulates cross-contamination and one that
assigns an abnormal amount of paralogs for the provided taxa.

**1. Replace duplicate sequences in receivers**

This mode is meant to simulate cross-contamination by replacing any duplicate
sequence with a randomly chosen sequence from another species.

**2. Assign a fixed number of sequences to receivers (abnormal paralogy frequency)**

## Dependencies

Python 3.5+ is assumed. Type `python3 --version` to see which version of Python
is on your system.

## Usage

Input alignments must be in the Fasta format. Put all of the alignments which
you want to run this script on into a single directory and provide the path to
that directory by typing the following into your command line.

```bash
./seqexchange --dir <alignment directory>
```

The script will automatically recognize these alignments, if their filetype
extension is in one of the valid formats (see below).

| Filetype | Extensions recognized by SeqExchange |
| -- | -- |
| FASTA files | `.fa`, `.fas`, `.fasta`, `.fna`, `.faa`, `.fsa`, `.ffn`, and `.frn` |

You also need to provide a list of taxa which will receive sequences from
others. Receiver taxa are entered by using the `--receivers` flag:

```bash
./seqexchange --dir <alignment directory> --receivers Drosophila Caenorhabditis
```
