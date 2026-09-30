# Integrator: `learnings-review` v2 — REJECTED, NARROWLY. **Two tokens and four lines. Everything else is right.**

**Gate GREEN** at `049550f`: 145 launcher + 249 reduction, both `EXIT=0`, zero failures, DONE
marker; `pixi.lock` byte-identical to `exp`; ruff clean; clean fast-forward.

**v2 did what it said.** Its title — *"the review was right on every checkable point"* — is accurate,
and design confirms **five of six** corrections landed **accurately**, not merely present. The
second-order risk I briefed for (a fresh overclaim while correcting under a rejection) **did not
materialise**, and the v7 staleness I predicted **is not there** — section G already reads *"the whole
layer-(b) pre-run override was **removed** by the amendment-20 decompose (Slug A, `157821f`)"*, which
it had when written. **My prediction was wrong and design checked it fairly rather than adopting it.**

Three findings block. Each has a remedy already written and verified by the reviewer who found it:
**one token, one token, four lines.** Nothing else in this work order needs action.

---

## B1 — BLOCKING (one token). The B1 correction re-attached the wrong SHA

Section G item 3 now reads *"attempted at v6 **`e1c0d63`** and rejected at `d5733af`"*. I verified:
`e1c0d63`'s diff to `settings_editor.py` contains **zero** references to `_record_edit` or
`_session_edits` — it is *"items 2-5 — record corrections, sidecar pin, --timeout, LAYERS table"*
(18 insertions). The generation-3→4 attempt is **`3d71164`**, *"record a per-angle edit per CELL, not
per column"* (85 insertions, rewriting `_record_edit(self, name, row=None)` and moving
`_session_edits` to tuple keys).

**Why this is the correction's own failure.** My v1 remedy text was deliberately SHA-free
(*"item 3 reads 'attempted at v6, rejected at `d5733af`'"*). v2 changed the verb and appended the new
SHA while **re-attaching `e1c0d63` from the sentence it was replacing** (`"fixed at v6 e1c0d63"`). And
the new header asserts *"Every correction below is transcribed from `d5733af` and **re-verified
here**"* — but `d5733af` never names `e1c0d63`. A reader following the pointer to check "was
generation 3 attempted here?" finds items 2-5 and no recorder change: **a pointer into a real
artefact that does not say what the sentence claims**, in the one section the correction rebuilt.

**Fix: `e1c0d63` → `3d71164`.** Every *other* `e1c0d63` citation in the document (the D residue
rename, E's prose, F's timeout) is correct — those genuinely are items 2-5.

## B2 — BLOCKING (one token). The numpy half of the F1 fixture remedy is vacuous, and the comment says otherwise

`tests/unit/lr_reduction/test_save_config_json.py:63`, `:80-84`.

**`np.float64` is a subclass of Python `float`**, so `json.dumps` serialises it natively. Measured on
numpy 2.1.3:

```
np.float64  isinstance(float)=True   json.dumps -> OK
np.float32  isinstance(float)=False  json.dumps -> TypeError
np.int64    isinstance(float)=False  json.dumps -> TypeError
np.bool_    isinstance(float)=False  json.dumps -> TypeError
```

So `config.IncidentTheta = np.float64(4.0)` needs no conversion, and
`isinstance(payload["IncidentTheta"], float)` is a tautology of numpy's type hierarchy rather than a
property of the code. The in-test comment *"A PosixPath and a numpy scalar both had to be converted on
the way out; neither is JSON-serializable"* is **false for the numpy half**, as is the commit body's
*"asserts both survive to JSON"*.

**The consequence, which I reproduced myself:** the whole redness of "drop `make_json_safe`" rests on
the single `Path(` at `:60`. Drop `make_json_safe` **and** revert that one line to `str(...)` — the
exact drift the new comment warns about — and the suite is **9 passed**. The fixture is one line from
re-vacuity.

**Fix: `np.float64(4.0)` → `np.float32(4.0)`** (or `np.int64`/`np.bool_`). test-reviewer verified the
anchor then reds R7 *independently* of the `Path` line, which is what the fixture was for. `float32`
is also the better fidelity choice — `binary_processing.py:182` reads the PV through h5py, yielding
the file's dtype, and `float64` is the one numpy float that is JSON-native and therefore cannot pin
this branch. **Fix the comment in the same edit.**

### My own verification was insufficient, in a way worth naming

I ran "drop `make_json_safe` → **1 failed**" and reported the remedy as real. That was **true but
incomplete**: I verified the *aggregate* reds, never that *each anchor* reds independently. A
prescription claiming two anchors was satisfied by one, and aggregate redness proves only that **≥1**
anchor functions.

> **When a fix adds N anchors, mutate against each one independently. The aggregate reding is not
> evidence that any particular anchor works.**

This is amendment 16's one-row-per-independently-removable-clause rule applied to **test fixtures**,
and it is the same defect class I rejected T3 v7 for one cycle earlier — a guard that covers one
spelling while being believed to cover all. I held that to account and walked past this.

## B3 — BLOCKING (four lines). "Serialise before open" is pinned only for failures that happen *before* serialisation

`src/lr_reduction/new_reduction_from_file.py:496-500`; guard at
`tests/unit/lr_reduction/test_save_config_json.py:112-130`.

All three `bad` params raise **`AttributeError` at `:498` on `config.__dict__`** — the same line, the
same reason. **None reaches `json.dumps`.** So the parametrisation pins *"a non-config is rejected
before `open`"*, not the clause the implementation comment names (*"`json.dump` raising part-way
through left a truncated file behind"*). Side effect: the `TypeError` half of both `pytest.raises`
tuples is unreachable.

**Demonstrated:** mutation M1 keeps the reject-before-open gate and removes **only**
serialise-before-open → **9 passed**. With a real config carrying a value `make_json_safe` passes
through untouched (`np.bool_`, a `set`, a `datetime` — anything hitting the bare `else: return obj` at
`save_reduced_data.py:130`), the prior settings file went from **29 bytes of valid JSON to 546 bytes
of truncated output**, `prior settings intact: False`. Shipped v2: `intact: True`. **R6's harm,
reproduced by a mutation R6's guard does not catch.**

This is the same shape as v1's F2 — *a guard that pins a weaker, adjacent property than the one
claimed* — recurring **in the remedy for F2**. The assertion itself is right (`read_bytes() == before`
is byte-identity, exactly the claimed property); it is the **input set** that never reaches the
clause.

**Fix, verified:** add a fourth param that is a real config carrying a pass-through value (a four-line
module-level helper; `NRReductionConfig()` does no I/O at init). At HEAD: **10 passed**. Under M1:
**1 failed** on `[bad3]`, at the `read_bytes()` assertion. The existing `raises` tuple already admits
the `TypeError`.

---

## WHAT PASSED — verified, do not re-examine

- **Both v1 blockers are closed, by measurement.** I re-ran both previously-surviving mutations:
  dropping `make_json_safe` → **1 failed** (was 6-passed-survived); open-before-serialize + `unlink` →
  **3 failed** (was survived). Restore sha- and symbol-verified both times.
- **v2 is strictly additive** — the only deletion in the test file is the one fixture line that was
  upgraded. Every v1 function, assertion and param survives verbatim.
- **B5 is genuinely closed.** Both `example_nr_reduction.py` sites (`:220`, `:285`) now call
  `save_config_json`; the byte-identical `open()`-then-`json.dump` blocks are gone and `json`/`save_fn`
  are no longer imported there. A tree-wide search finds **no third copy** of that shape.
- **Five of six synthesis corrections are accurate**, checked against the artefacts: B1′'s description
  (the collapse at `settings_document.py:269`, the `elif` at `:534`, `LambdaMinUse = None` reaching
  `:452`), B2's masking direction, **B3's attribution** (`harness-hardening-learning.md §2` at
  `c6628ea` 2026-08-19 states it verbatim; `settings-editor-learning.md §8` appears 21 days later, one
  day *after* amendment 16's adoption, and is the diagnostic corollary), **B4's attribution** (the name
  is the plan's at `:953`, hedged "e.g."; `grep fifo_sidecar 8c7dfbb:todo.md` returns nothing), and F's
  timeout. The five folded advisories all landed.
- **My v1 `__dict__` statements were superseded correctly** — `Spath` *does* have a setter, so v1 was
  wrong; the real trap is `base_path`, and the re-derive-vs-freeze reasoning holds.
- **No v7 staleness**, and the `todo-untested-new-workflow-launcher-path` annotation landed with the
  corrected B6 framing.

## ADVISORY

1. **"the count is 1" is true only of the open-then-dump family.** `SettingsDocument.save`
   (`settings_document.py:305`) serialises the *identical* payload — `to_dict()` is literally
   `dict(self._config.__dict__)` — but writes it atomically (`mkstemp` + `fsync` + `os.replace`) and
   refuses symlinks. So the behaviour has two implementations and **the surviving public one on the
   scientist path is the weaker**: `open(path, "w")` still truncates at open, so ENOSPC or SIGKILL
   between open and write still destroys the prior file. Pre-existing, not worsened — but one sentence
   in the code section would stop "the count is 1" reading as closure of the class that section is about.
2. **G item 4 drops a condition both its sources carry.** The collapse fires only
   `if field.optional_list and all(entry is None …)`; `d5733af` and §21 both say *"when the last
   populated cell is cleared"*. As written it implies one gesture suffices on a multi-angle column,
   which the same sentence's "two-gesture variant" contradicts.
3. **"4 of 5 red" understates its own frame** — the bullets above it disclose two further mutations at
   unenumerated sites, both survivors. Post-`0f35b83` the honest score is **6 of 7 red, one disclosed
   survivor**. Amendment 16's frame clause applied to this review's own ledger.
4. **The two new `example_nr_reduction.py` delegations are unpinned by construction** — reverting both
   to inline `json.dump` gives 9 passed, and nothing outside the file references it. test-reviewer
   ruled Advisory for consistency with its own v1 F5 on the structurally identical `:152` survivor: no
   cheap *behavioural* guard exists for either, and the cheap alternatives are the adjacent-property
   anti-pattern this review rejects. **I agree.** Extend the todo annotation to name these two sites.
5. **Fixture fidelity is right for one field, stale for three** — `file_batch.py:615-620` puts `Path`
   into `DBpath`, `Spath` **and** `NEXUSpathRB`, while the fixture uses `str()` for three of four.
   Not a gap today; it is the drift surface B2 makes load-bearing.
6. **Header mis-assigns two of four** — B5/B6 came from `2136a8f`'s **own** body, and B3 was an
   independent attribution error, not inheritance. The stated root cause (a self-report is a claim,
   not a record) covers both bodies and is right; only the commit named is too narrow.
7. **One event behind, fairly.** My v7 rejection `90fd351` (14:17, after the synthesis) adds a fresh
   instance to **this document's own class A** — the layer-(b) pin is blind to `dataclasses.replace`
   (203 green under the mutation) — and to class E, its false docstring. It also records that
   `plans/settings-ui-override-plan.md` has **zero** occurrences of `ipts`/`forget`/`experiment`, so
   Slug B's plan owes the experiment-bound science invariant as well as cell-level layer authority.
   Fold at the next revision; nothing currently written is wrong.
8. Record precision: the body says `:68` stores `datapath` "straight into `_NEXUSpathRB_override`" —
   it goes through the public setter, which happens to be a bare store. Non-material.
