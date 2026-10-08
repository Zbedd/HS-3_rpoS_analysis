# RpoS alignments

Run `python scripts/reproduce.py --only alignments` from the project root to
stage the frozen inputs and generate the four manuscript alignment panels,
diagnostic tables, findings, and run manifest. Outputs are written under
`outputs/rpos_sequence/rq2/rpoS/` in the selected reproduction workdir.

For a prepared local cache, use `python -m seq.rpos_alignments --workdir DIR`.
Additional domain cartoons and sequence logos require `--exploratory`.
See [methods](methods.md) for the retained comparisons and citations.
