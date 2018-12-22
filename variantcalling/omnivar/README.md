# SwitchVar

Switchvar is a local assembly-based variant caller for short-read DNA. It is able to detect long indels and deletions that are missed by alignment-based callers, while performing similarly to alignment-based callers on short indels and SNPs. Switchvar performs local reassmbly on strategically selected small intervals of the reference genome using a read-aware DeBruijn graph. Its switchyard algorithm determines possible haplotypes based on read continuity through paths in the Debruijn graph. Each haplotype is scored based on its read coverage characteristics, and a genotype is called based on the scores.

### Installing

From the directory containing setup.py, run "pip install ."

## Getting Started

More to come.

## Running the tests

Running the script "runtests.sh" will run all tests found below your current working directory.

## Versioning

We use [SemVer](http://semver.org/) for versioning.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
