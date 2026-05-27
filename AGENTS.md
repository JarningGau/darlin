# Repository Guidelines

Contributor guide for **darlin**, a Python 3.11 project managed with [Pixi](https://pixi.sh/). The tree is still small; conventions below describe the current setup and expected layout as code grows.

## Project Structure & Module Organization

- **Root:** `pixi.toml` and `pixi.lock` define the workspace, channels (conda-forge, bioconda), and `linux-64` platform.
- **Environment:** Pixi installs dependencies into `.pixi/` (ignored by Git except `.pixi/config.toml` if present).
- **Data:** `data/` and `old/` are listed in `.gitignore`—keep large or sensitive datasets out of version control unless you intentionally add exceptions.
- **Source & tests:** Put application code under a clear package name (e.g. `darlin/` or `src/darlin/`). Tests live under `tests/` with parallel module paths to the code they cover; larger or shared fixtures may sit under `tests/data/` (or similar) instead of the repo-root `data/` tree.
- **Docs:** Keep command references in `docs/` aligned with CLI behavior, especially when output filenames or tabular output schemas change. For user-visible behavior or output changes, add a short entry under **Unreleased** in `docs/CHANGELOG.md`.

## Build, Test, and Development Commands

| Command | Purpose |
|--------|---------|
| `pixi install` | Sync the locked environment from `pixi.lock`. |
| `pixi run <task>` | Run a task defined under `[tasks]` in `pixi.toml`. |

Common tasks in this repo:

| Task | Command | Purpose |
|------|---------|---------|
| smoke | `pixi run smoke` | CLI help/argparse smoke checks |
| compile | `pixi run compile` | `compileall` over `src/` and `tests/` |
| test | `pixi run test` | run pytest |

After adding `[tasks]` entries, document them in this file or in `README.md` so others can discover them quickly.

## Coding Style & Naming Conventions

- **Language:** Python **3.11** (see `[dependencies]` in `pixi.toml`).
- **Indentation:** Four spaces; no tabs for Python.
- **Naming:** `snake_case` for functions and modules; `PascalCase` for classes; constants in `UPPER_SNAKE_CASE` where appropriate.
- **Formatting & linting:** Not pinned in-repo yet. When you add tools (e.g. Ruff, Black), list them in `pixi.toml` and prefer `pixi run` tasks so everyone uses the same versions.

## Testing Guidelines

- **Framework:** Prefer **pytest**; name files `test_*.py` or `*_test.py` and keep them under `tests/` (or colocated next to modules if you adopt that style for a specific area).
- **Running tests:** Use `pixi run test` once that task exists, or `pytest` from an activated Pixi environment.
- **Coverage:** No global requirement yet; add a sensible threshold when the codebase justifies it.

## Commit & Pull Request Guidelines

- **Commits:** Until there is shared history on the default branch, treat this document and `pixi.toml` as the source of truth for conventions. **Commit messages must follow [Conventional Commits](https://www.conventionalcommits.org/):** a type prefix (`feat`, `fix`, `chore`, `docs`, `refactor`, `test`, `ci`, etc.), optional scope in parentheses, then a short imperative description (e.g. `feat(cli): add dry-run flag`). Use the body for context when the subject alone is not enough.
- **Pull requests:** Describe what changed and why; link tracking issues when applicable. For user-visible or workflow changes, note how to verify (commands run). This repo is currently CLI/library-oriented.

## Lockfile & Merges

`pixi.lock` is configured in `.gitattributes` for merge safety. Resolve lockfile conflicts by regenerating (`pixi lock` / `pixi install`) rather than hand-editing when possible.

## Version bumps

The published package version is `__version__` in `src/darlin/__init__.py` (wired into `pyproject.toml` for editable installs and builds). Documented releases and commit anchors live in `docs/CHANGELOG.md`.

When cutting a new release:

1. **Changelog:** Add `## [x.y.z] - YYYY-MM-DD` above **Unreleased** in `docs/CHANGELOG.md`, move finished bullets out of **Unreleased** into that section (or rewrite for accuracy), then reset **Unreleased** for ongoing work.
2. **Version constant:** Bump `__version__` in `src/darlin/__init__.py` to `x.y.z` so it matches the new changelog heading.
3. **Git tag (optional):** After the release commit is on the branch you ship from, create an annotated tag `vx.y.z` if your workflow uses tags.

Use [SemVer](https://semver.org/) intent: during `0.x`, treat incompatible CLI or output-schema changes as at least a **minor** bump unless you explicitly document a patch-only policy; call out breaking changes in the changelog.

## Agent-Specific Instructions

- Read `pixi.toml` before suggesting dependencies or Python version changes.
- Do not commit contents of `data/` or `.pixi/` unless the user explicitly wants them tracked.
- Prefer adding Pixi **tasks** for repeatable commands instead of one-off shell instructions in docs alone.
- When changing scRNA pipeline output schemas, update `tests/test_cli_scrna.py` and the relevant command docs together so the documented outputs stay synchronized with the CLI.
