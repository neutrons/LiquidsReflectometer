# Release Process for the `exp` Branches

Abridged from the NDPD [Software Release Process](https://ornl-neutrons.atlassian.net/wiki/spaces/NDPD/pages/23363585/Software+Release+Process),
specialized for the experimental line of LiquidsReflectometer. The `exp` line follows the same
three-branch model as `next`/`qa`/`main` but with its own branches, tags, and conda labels, so it
never collides with the production line.

## Branches

| exp line   | mirrors | role               | conda label | publishes on        | deploy                                             |
|------------|---------|--------------------|-------------|---------------------|----------------------------------------------------|
| `exp-next` | `next`  | development        | `exp-next`  | every push          | automatically on push                              |
| `exp-qa`   | `qa`    | release candidates | `exp-qa`    | `exp*rc*` tags only | manually update tag and run pipeline in GitLab     |
| `exp`      | `main`  | releases           | `exp`       | `exp*` tags only    | manually update tag and run pipeline in GitLab     |

- Merges flow `exp-next` → `exp-qa` → `exp` and are fast-forward only (`--ff-only`).
- Every push to `exp-next` also triggers the GitLab deploy pipeline for `lr_reduction_dev-exp`. The workflow does not deploy `exp-qa` or `exp` automatically.

## Versions and tags

Tags are named `exp<version>`; versioningit strips the `exp` prefix to get the package version.

| tag            | package version | conda label |
|----------------|-----------------|-------------|
| `exp3.0.0rc1`  | `3.0.0rc1`      | `exp-qa`    |
| `exp3.0.0`     | `3.0.0`         | `exp`       |

- Use the full three-component form, `exp<major>.<minor>.<patch>` and `exp<major>.<minor>.<patch>rc<n>`.
- The exp line starts at 3.x so its versions sort above the `2.10.0.dev` builds already on the `exp` label.
- Between tags, `exp-next` builds as `<next minor>.dev<timestamp>+g<sha>`, for example `3.1.0.dev20261026202745+gb712a3f` after `exp3.0.0rc1`.
- Never push `v*` tags from an exp branch. Those belong to the production line.

The publish job derives the conda label from the tag name (`rc` → `exp-qa`, otherwise `exp`) and
refuses to publish unless the tagged commit is already on that branch at `origin`. Push the branch
before the tag.

## Marking a release candidate

### 1. Sync `exp-next` with the remote

```bash
git fetch -p -P -t --all
git checkout exp-next
git rebase -v origin/exp-next
```

### 2. Make sure `exp-qa` and `exp` are already in `exp-next`

```bash
git merge --no-edit origin/exp-qa
git merge --no-edit origin/exp
git diff origin/exp-next   # expect no output
```

Both merges should answer `Already up to date`. If they did not, push `exp-next`.

### 3. Fast-forward `exp-qa` to `exp-next` and push

```bash
git checkout exp-qa
git rebase -v origin/exp-qa
git merge --ff-only origin/exp-next
git push origin exp-qa
```

### 4. Tag and create the pre-release

Either let the GitHub release create the tag:

```bash
gh release create exp3.0.0rc1 --title "exp3.0.0rc1" \
  --target exp-qa --generate-notes --prerelease
```

or tag manually and push the tag:

```bash
git tag exp3.0.0rc1
git push origin exp3.0.0rc1
```

In the GitHub UI: Releases → Draft a new release → tag `exp3.0.0rc1`, target `exp-qa`,
Generate release notes, check **Set as a pre-release**, Publish.

The tag push runs the workflow, builds `lr_reduction 3.0.0rc1`, and publishes it to Anaconda under the `exp-qa` label.

## Marking a release

Only promote a tagged release candidate. Changes present in `exp-next` but not in `exp-qa` go through a release candidate first.

### 1. Sync `exp` with the remote

```bash
git fetch -p -P -t --all
git checkout exp
git rebase -v origin/exp
```

### 2. Fast-forward `exp` to the release candidate and push

```bash
git merge --ff-only exp3.0.0rc1
git push origin exp
```

### 3. Tag and create the release

```bash
gh release create exp3.0.0 --title "exp3.0.0" \
  --target exp --generate-notes --latest=false
```

or manually:

```bash
git tag exp3.0.0
git push origin exp3.0.0
```

In the GitHub UI, target `exp` and leave **Set as the latest release** unchecked so the production
line's release stays the repository's "Latest".

The tag push builds `lr_reduction 3.0.0` and publishes it under the `exp` label.

## Anaconda labels

- `exp-next`: development builds from every push. The workflow keeps the five most recent and removes older ones.
- `exp-qa`: release candidates. Builds accumulate; remove old ones manually if needed.
- `exp`: releases.

The deploy environments use packages from their respective label, e.g. `lr_reduction_exp_qa` gets `lr_reduction` packages with the label `exp-qa`.

## Configuration (`pyproject.toml`)

```toml
[tool.versioningit.vcs]
method = "git"
default-tag = "0.0.1"
match = ["exp[0-9]*"]      # only exp<version> tags; v* and exp-review-* tags are ignored

[tool.versioningit.tag2version]
rmprefix = "exp"           # exp3.0.0rc1 -> 3.0.0rc1

[tool.versioningit.next-version]
method = "minor"
```
