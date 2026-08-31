# Integrator review gate — `scaling-factor-path-anchor` v2 REJECTED

Supersedes the v1 rejection at `aac75a9` (`git log -- todo.md`).

**Disposition: blocking, per the plan's own declaration.** `test-reviewer` is
declared blocking for this bug-fix phase and returned one blocking finding. No
Integrator departure was needed. Retry budget: attempt 2 of N = 3.

**The fix is two lines.** Everything else below is advisory.

## v2 is very good work — this is the larger part of the truth

All three v1 blocking findings are closed, and every one was verified
**red-on-revert** by mutation rather than accepted on report:

| reverted | result |
|---|---|
| the `_resolve_scaling_factor_file` call in `read_template` | RED |
| `raise FileNotFoundError` → `print` + `return workspace` | RED |
| reverse the candidate precedence (anchored first) | RED |
| delete the as-given candidate | RED |
| absolutise the fixture, anchoring left intact | RED (the new `isabs` precondition) |
| unresolvable → return the anchored candidate | RED |

**All six advisories were taken**, including two I would have accepted without:
the upward directory walk is **gone entirely** (one anchored candidate, magic
depth constant removed) with the reasoning in the comment; `logger.notice` on
substitution and `logger.warning` listing candidates on failure; the no-match
unscaled path logged (*"the reflectivity comes out plausible and wrong"*);
`meta_data["scaling_factor_file"]` recorded; and the forked `read_template` in
`new_reduction_from_template.py` anchored identically so the two readers cannot
diverge.

**And the change did something the commit body undersells.** Because the six
fixtures now declare `sf_197912_Si_auto.cfg` (template-relative) rather than
`data/sf_...`, the anchoring is load-bearing **under `cd tests/` as well**.
Measured: with the anchoring deleted and cwd `tests/`, `read_template` returns
an unresolvable path. So deleting the fix now turns **the campaign gate itself**
red — not merely bare pytest. v1's guard worked only because its tests
deliberately `chdir`'d away; v2 makes it unconditional. **That retires my gate's
structural blindness to this class**, which is worth stating plainly in the PR.

**The cwd leak was masking 11 tests**, all now fixed — measured on base
(9 real failures, reproducing the human's report exactly), base + the `chdir`
fix alone (15), and HEAD (**0**). All are discoveries, not regressions, and the
commit body's arithmetic checks out.

## B4 — BLOCKING: `tests/test_time_resolved.py` is not green standalone

The plan's v2 acceptance criteria require the two spot-check files *each* green
from the root. Measured on this clone, repo root:

| invocation | result |
|---|---|
| `pytest tests/test_scaling_factors_workflow.py tests/test_time_resolved.py` (the commit's spot-check, verbatim) | **7 passed** |
| `pytest tests/test_scaling_factors_workflow.py` alone | **5 passed** |
| `pytest tests/test_time_resolved.py` **alone** | **2 failed** |

```
ValueError: Invalid value for property Filename ... "REF_L_198413" not found
src/lr_reduction/time_resolved.py:120  api.LoadEventNexus("REF_L_%s" % meas_run_30Hz)
```

**Mechanism, confirmed by the Integrator.** `tests/test_time_resolved.py` is the
**only** Mantid-loading module in the suite that never sets the facility:

```
tests/test_time_resolved.py             0 occurrences of default.facility
tests/test_reduction.py                 1
tests/test_dead_time.py                 1
tests/test_scaling_factors_workflow.py  1
```

and `amend_config` (`src/lr_reduction/utils.py:63-69`) sets `default.facility`
**only when `new_config` is passed** — these tests pass `data_dir` only. So the
module is green only when co-collected with a module that sets the facility.

**This is the same species as B1, one cycle later**: a green produced by another
test's process-global side effect, reported as though it were a property of the
file. Not a regression — it was red standalone on base too, for a different
reason — but the v2 claim that it is now green standalone is wrong, and
eliminating exactly this accident is what the retry existed to do. Letting the
shape through would make the gate decorative.

**Fix (verified by the reviewer: 2 passed alone from the root):**

```python
import mantid.simpleapi as mtd_api
mtd_api.config["default.facility"] = "SNS"
mtd_api.config["default.instrument"] = "REF_L"
```

Better: hoist those into `tests/conftest.py` as a session-scoped autouse fixture
so no module can be silently coupled again. Then correct the commit/PR body —
the spot-check line must read as **two separate invocations**, with the
standalone number stated.

## Strongly recommended in the same commit — closes the class, not the instance

1. **Delete the now-vestigial `chdir`** at `tests/test_reduction.py:235`. Every
   path in `test_reduce_functional_bck` is already absolute; the reviewer
   deleted the line and ran the test alone from the root — **1 passed**.
   Converting the leak to `monkeypatch.chdir` mitigates one instance; deleting
   it removes the class, and the suite then contains zero cwd mutations.
2. **Add a cwd-leak detector** — nothing currently catches a reintroduced leak
   (reverting `monkeypatch.chdir` → bare `os.chdir` leaves the suite green,
   because v2 removed every cwd dependence and thereby made the leak inert):

   ```python
   @pytest.fixture(autouse=True)
   def _no_cwd_leak():
       before = os.getcwd()
       yield
       assert os.getcwd() == before, f"test leaked cwd: {before} -> {os.getcwd()}"
   ```

   Four lines converts B1's whole failure class from "found by a human months
   later" into "found by the next test in the file".
3. `tests/test_template_paths.py:1-12` — the guard file's own docstring is
   already stale: it quotes the pre-v2 fixture path, says the gate "structurally
   cannot see" the defect (which v2's fixture change made false), and says "Both
   tests below" where there are now four.
4. `:39` — the precondition uses `re.search`, which reads the **first**
   `<scaling_factor_file>`, while the test reads sequence 7. Use `re.findall`
   and assert none is absolute.
5. `:97`, `:109` test the private `_resolve_scaling_factor_file` directly and
   build XML that is never parsed. Route both through `read_template` — same
   test, pins the public seam.

## For the PR body / separate slugs

- **A second live global-state leak, same species, load-bearing on 14 tests.**
  `amend_config` (`utils.py:86-99`) takes a backup of `datasearch.directories`
  but never appends the key to `modified_keys`, so the `finally` never restores
  it — the backup is dead code. The reviewer added the missing line and got
  **14 real failures**, because tests call `process_from_template_ws` *outside*
  their `with amend_config(...)` block and rely on the leak leaving the nexus
  dir on Mantid's global search path. Out of scope here, but it means the PR
  should claim **cwd**-independence, not independence generally. Deserves its
  own slug with the seriousness B1 got.
- **Seven tracked templates, not six.**
  `tests/data/template_with_instrument_settings.xml` still declares
  `test/data/sf_186529_Si_auto.cfg` — note `test/`, resolvable from nowhere.
  Inert today, but it now emits the new "Could not resolve" warning on every
  read, and it is the last template on the old convention.
- `meta_data["scaling_factor_file"]` is set unconditionally, including on the
  "Skipping scaling factor" branch where nothing was read.
- The `new_reduction_from_template.py` parity line is untested (that module has
  **zero** test references), and unlike `template.py`'s writer its
  `write_template` serialises the resolved absolute path into
  `REFL_<seq>_template_new.xml`, giving a saved template with mixed path
  conventions.
- `src/lr_reduction/template.py:466` has a cwd-relative write in **product**
  code, currently dead behind `OUTPUT_NORM_DATA = False`. One line, worth doing
  while the context is loaded.
- `time_resolved.py:257,374` bare `except:` still swallow the new
  `FileNotFoundError` — correctly flagged in the commit body, and confirmed not
  to create a vacuous green (both tests index `reduced[0]`, so an all-slices-
  failed run raises `IndexError` rather than passing on empty data).

## Gate result — both invocations green

- **Charter** (`pixi run test-reduction`, cwd `tests/`): `test-launcher` **10
  passed**; reduction **111 passed**, exit 0, zero failures.
- **Repo root** (`pixi run pytest`): **111 passed**, exit 0, zero failures.
- **And green for the right reason.** I re-applied the discriminator that
  exposed v1's fake green: `test_scaling_factors_workflow.py` alone from the
  root and the five previously-failing tests by node id — **5 passed and 5
  passed**, where v1 gave 5 failed and 5 failed. Selection-independence is
  restored on the cwd axis. B4 is the one file where a *different* global
  coupling remains.

Note for the record: the reviewer disclosed that its probes overwrote
`/tmp/sf_197912_Si_test*.cfg` between ~22:50 and ~23:15 while exercising those
five tests. My gates passed, so nothing here is attributable to it — but it is
an independent reproduction of the shared-`/tmp` race, and the write side of
those five functions (`output_dir = "/tmp"`, four sharing one filename) is still
unfixed while v2 fixed their read side. That also violates the repo's
Secure Temporary Files rule (predictable, world-readable path).

## Suggested v3 order

1. **B4** — the two `mtd_api.config` lines, preferably as a conftest autouse
   fixture, and correct the spot-check claim in the commit/PR body.
2. Delete the vestigial `chdir` and add `_no_cwd_leak` — together they close the
   class B1 belonged to.
3. The `test_template_paths.py` docstring, `re.findall` precondition, and
   routing the two helper tests through `read_template`.

## Addendum (reviewer's completed matrix) — one coverage gap worth closing in v3

The blocking finding above is unchanged. Two items the fuller report adds:

- **Nothing pins the symptom the human actually reported.**
  `test_scaling_factor_missing_file_raises` pins the *leaf function*. What was
  observed in the field was the failure surfacing **through the pipeline** —
  `TypeError: cannot unpack non-iterable EventWorkspace object` four frames
  above, at `template.py:385`. A `pytest.raises(FileNotFoundError)` around
  `process_from_template_ws` with an already-loaded fixture workspace would
  close the loop at the level the bug was seen, and would fail if a future
  refactor reintroduced an unpack-mismatch anywhere between the two. Cheap, and
  it is the one test that maps 1:1 onto the report that started this slug.
- `tests/test_reduction.py:90, 125, 179` — the three `cleanup_partial_files`
  calls are dead (`process_from_template_ws` writes no partial files), and three
  tests gained a `tmp_path` parameter solely to feed them. Deleting the calls and
  the parameters is the cleaner end state. (Unrelated and genuinely good:
  repointing the real `output_dir`s at `tmp_path` is a real isolation win —
  those tests previously wrote partials and `REF_L_*_auto_template.xml` into
  `tests/data/`, which is why working clones accumulate the 16 gitignored
  byproducts that make them diverge from a fresh checkout.)

Also recorded: the reviewer's probes overwrote the shared
`/tmp/sf_197912_Si_test*.cfg` at 23:12–23:13. Both of my gates passed, so
nothing in the result above is attributable to that — but it is an independent
reproduction of the shared-`/tmp` race and strengthens the case for folding the
**write** side of `test_scaling_factors_workflow.py` into the parked
`test-tmp-isolation` slug.
