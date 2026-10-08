# `src/` — the importable code

Everything here is installed as a package, once per environment:

    .venv/bin/pip install -e .

That is what makes `paths`, `viz`, `assays`, `seq` and `figures` importable
from any working directory. Nothing in this repository manipulates
`sys.path`, and nothing counts `..` from its own file.

```
paths.py           REPO_ROOT, data(assay), outputs(assay)
viz/
  style.py         colours, type, marks, chrome, export
  layout.py        panel arrangement for composed figures
assays/<assay>/    config, io, stats, panels, figure
seq/               sequence-validation workstreams + the shared seqfetch cache
figures/           manuscript figures composed from assay panels
```

CLI drivers live in `../scripts/` and hold argparse, reporting and the
orchestration that decides which figures an assay emits. A driver composes;
it does not compute or draw.

`scripts/reproduce.py` stages the frozen manuscript inputs in an isolated run
directory and calls the existing drivers. `HS3_WORKDIR` redirects their data
and output paths for that run; `paths.SOURCE_ROOT` still locates this checkout.
See the root README for the public entry point.

## The assay contract

Every `assays/<assay>` package exposes the same four things:

| module   | holds                                                              |
|----------|--------------------------------------------------------------------|
| `io`     | `load…()` returning a tidy frame — the only place units, dilutions and censoring are defined |
| `stats`  | anything numeric, importable without selecting a plotting backend  |
| `panels` | functions that draw one panel into an axes the caller owns, creating no figure and writing no file |
| `figure` | the assay's own standalone figures, each composing its panels      |

The panel function is the load-bearing part: it is what lets one piece of
code serve both the assay's own figure and a manuscript figure, so the two
cannot disagree about what a panel looks like.

Where a panel and a composed figure have to agree on a layout, one of them
owns it and hands it over rather than both deriving it. `assays/cr/panels.py`
is the worked example — `level_positions` lays out the levels columns once,
and the levels panel, the group labels beneath it and the colony strip above
it all read from that one result.

## Dependency direction

`figures` imports `assays`. `assays` never imports `figures`. A figure that
spans two assays therefore has a home that neither assay owns.

## Strain naming

The names in the data files are not the names in the paper, and they are not
meant to be. Everything under `data/` uses the names the strains were built
and recorded under; `viz/style.py` translates those into what a reader sees,
and that translation is the only place the paper's naming exists.

| in `data/` and the configs | in every figure and the paper | what it is |
|---|---|---|
| `WT` | **Parental** | the original unmarked isolate, carrying no cassette |
| `kanR` | **WT** | the kanamycin-marked strain, and the reference every comparison is made against |
| `rpoS`, `rpoS::ΩkanR` | ***rpoS*⁻** | the *rpoS* mutant |

"WT" in an HS-3 figure denotes the kanamycin-marked control. The manuscript
describes a random insertion in this control and a separate, targeted insertion
in the *rpoS* mutant. The unmarked isolate is labelled Parental and is recorded
separately in the source data.

`fs.lesion(gene)` builds the superscript-minus form; PA14 is its own
organism's wild type, so it too is labelled WT.

To rename a strain, edit `LABELS` in `viz/style.py` and nothing else. The
raw names stay put, so colours, markers and every CSV are untouched by a
renaming — only the text on the figures moves.

## Figure style

`viz` is the single source of truth for how a figure looks and how its panels
are arranged. Every plotting module does `import viz as fs` and calls
`fs.use()` once from its driver; nothing sets rcParams, colours or figure
sizes anywhere else. The categorical palette is validated for colour-vision
deficiency — the derivation and the re-validation command are in the
`viz/style.py` docstring.

`viz.layout` covers the multi-panel case: `stack()` builds a page from one
dict per lettered panel, each mapping row names to heights in **inches** (so a
row keeps its printed height when a neighbour changes) and taking the gaps
between them from `style`; `split()` subdivides a row into side-by-side
panels; `tag_panels()` letters them, and freezes the layout as it does, so it
goes last; `collect_legend()` merges a figure-wide key. Composed figures use
constrained layout; `tight_layout` does not solve the nested gridspecs they
need.

## Adding to the tree

A new package needs an `__init__.py` or setuptools will not find it. A
non-`.py` file read at runtime needs a `package-data` entry in
`pyproject.toml`, and a new top-level module needs a `py-modules` entry.
An editable install hides both omissions until someone builds a wheel:

    .venv/bin/python -m build --wheel --outdir /tmp/w .

## Changing a figure

Every pipeline is byte-reproducible under `SOURCE_DATE_EPOCH`, so a change
meant to alter nothing can be shown to alter nothing. Delete the outputs
first — a stale file that was never rewritten passes a hash check it never
earned.

    export SOURCE_DATE_EPOCH=1700000000
    rm -f ../outputs/<assay>/*.png ../outputs/<assay>/*.pdf
    .venv/bin/python ../scripts/<driver>.py
    # re-hash and diff against the hashes taken before the change

State up front which figures a change *should* alter, so the diff is a list
you predicted rather than one you discovered.
