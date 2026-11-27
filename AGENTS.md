# bacass: agents

nf-core/bacass is a Nextflow (DSL2) bioinformatics pipeline for simple bacterial assembly
and annotation. It assembles short reads, long reads, or a mixture of both (hybrid
assembly).

For every session, fetch and follow the instructions at
https://raw.githubusercontent.com/nf-core/agents/main/resources/pipeline/AGENTS.md before
editing anything. That file governs code style, prose, module and subworkflow structure,
testing, git flow and commit practice, and it wins any conflict with this file. Repository
conventions are described in `docs/CONTRIBUTING.md`.

## Commands

Requires Nextflow, nf-core tools, and nf-test.

Run pipeline tests (small test profiles, requires Docker/Singularity/Conda):

```bash
nf-test test --tag pipeline --profile +docker --verbose
```

Run a single test file:

```bash
nf-test test tests/default.nf.test --profile +docker --verbose
```

Regenerate snapshots:

```bash
nf-test test --tag pipeline --profile +docker --verbose --update-snapshot
```

Lint with the nf-core community rules:

```bash
nf-core pipelines lint .
```

Rebuild `nextflow_schema.json` after changing params in `nextflow.config`:

```bash
nf-core pipelines schema build
```

Bump the minimum required Nextflow version:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

Pre-commit hooks (prettier, trailing whitespace, nextflow-lint) are declared in
`.pre-commit-config.yaml` and run in CI in `.github/workflows/linting.yml`.
