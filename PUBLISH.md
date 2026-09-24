# Publishing txmod 0.1.0

Everything in this tree is ready. A few things are deliberately left to you,
because they carry your identity or an account you own.

## 0. Fill in the placeholders first

| File | Line | What to put |
|---|---|---|
| `pyproject.toml` | `authors = [{ name = "TxMod developers" }]` | the real author list |
| `CITATION.cff` | `authors:` and `repository-code:` | the same author list, and the repository URL from step 2 |
| `README.md` | `## Citation` -> `[PLACEHOLDER: citation on acceptance]` | the manuscript citation, or a preprint DOI if you post one |
| `README.md` | `<repo-url>` (Install section, four places) | the repository URL from step 2, e.g. `https://github.com/YOUR_USER/txmod.git` |
| `paper/README.md` | `Zenodo (DOI: [PLACEHOLDER])` | the DOI of the processed-tables deposit (the data deposit, not this repository's archive) |

## 1. Make the first commit under your own identity

    cd txmod_pkg
    git init -b main
    git add -A
    git -c user.name="YOUR NAME" -c user.email="YOUR EMAIL" \
        commit -m "txmod 0.1.0: transcript-resolved 3' UTR mutation interpretation"

Or set them globally first (`git config --global user.name ...`) and drop the `-c` flags.

## 2. Create the repository and push

With the GitHub CLI:

    gh repo create txmod --public --source=. --remote=origin --push

Or by hand: create an empty public repository named `txmod` on github.com, then

    git remote add origin https://github.com/YOUR_USER/txmod.git
    git push -u origin main

The CI workflow (`.github/workflows/tests.yml`) runs on the first push: pytest on
Python 3.9 / 3.11 / 3.12, plus a CLI smoke test on the bundled example files.

## 3. Tag the release, then archive it on Zenodo

    git tag -a v0.1.0 -m "txmod 0.1.0"
    git push origin v0.1.0

Then on Zenodo: log in with GitHub, switch the `txmod` repository on under
Settings -> GitHub, and publish a GitHub **release** from tag `v0.1.0`. Zenodo
picks that release up automatically and mints a DOI. `CITATION.cff` is what it
reads for author and title metadata, which is why step 0 matters.

## 4. What to send back

Two strings finish the manuscript:

- the repository URL -> Code Availability, `[PLACEHOLDER: repository URL]`
- the Zenodo DOI -> Code Availability, `archived at Zenodo (DOI: [PLACEHOLDER])`

A third Zenodo DOI is still needed for the processed intermediate tables
(Data Availability); that is a separate upload, not this repository.

## State of the tree as bundled

| | |
|---|---|
| tests | 46 passed, verified in a clean virtual environment |
| CLI smoke test | `txmod prepare` on the bundled example produced 6 variant-transcript pairs |
| version | 0.1.0 consistent across `pyproject.toml`, `src/txmod/__init__.py`, `CHANGELOG.md` |
| licence | MIT, `LICENSE` present |
| `paper/` | all manuscript figures re-render from the Zenodo tables and match the submitted files |
| earlier working names | none left in the tree (`txm6a`, `IsoPerturb`) |
