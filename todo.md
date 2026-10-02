# todo.md — Integrator rejection, `editor-load-fidelity` v1 @ 9fe4184 (review gate: design domain, harm clause)

**Verdict: REJECT — one blocking finding (demonstrated reachable harm).** Gate cycle: v1 → `review/editor-load-fidelity`
(retry 1 of N=3). Everything else passed; the defect is in the plan's prescription (A3 / §3 type-table row
"scalar `1`/`0` → loaded `True`/`False`, saved `true`/`false`"), so the fix starts with a plan revision.

## What passed (do not redo)

- **Gate** `pixi run test-reduction` from the subject root, analysis clone 2 (deploy-target class), 14:36–14:47 EDT:
  launcher **155 passed**, reduction **281 passed**, 699 warnings (= the base), **exit 0**. Matches the GREEN commit.
- **Diff scope**: the plan's five files only; no `plans/`, `todo.md` or battery.
- **Test reviewer: PASS.** Every §6 row present and type-safe (`True == 1` rule held line by line); M1, M2, M5, M9, M12,
  M22 re-executed independently on a scratch copy → exactly the recorded counts; others consistent by reading.
- **§8.5 deployment-shaped acceptance: PASS.** `ledger/scripts/editor-real-file-roundtrip.py` (sha256 `775c8c86f9cdb185`)
  → 7/7, exit 0, on six real reducer-written `*_settings.json` (IPTS-37957, -36119, -36574, -36681, -37096, -38511,
  `shared/reduced/`) and `IPTS-36119/…/new_reduction/REFL_234198_combined_autoreduction.dat`; 0/4 at the base (RED).
  Offscreen Load in the real tab (`ledger/scripts/editor-offscreen-load-probe.py`, sha256 `c993e19ec180cab8`): 3/3 at the
  tip (cells `true`, record read-only), FAIL at the base.

## BLOCKING — B-1: a hand-written `"useGravity": 1` loads silently as `True` and saves as `true`, turning gravity correction ON

**Consumer.** `src/lr_reduction/nr_reduction_calc.py:1079` — `if self.config.useGravity is True:` — an **identity** test.
The reducer's own loader does no type conversion (`new_reduction_from_file.py:439-449`, `json_to_config` = `setattr`),
so a file holding `1` reduces with gravity correction **off**.

**Reproduction** (Integrator, 9fe4184, analysis clone 2; temp files only):

```python
src.write_text(json.dumps({"useGravity": 1}))
json_to_config(json.loads(src.read_text())).useGravity is True        # False  -> reducer: gravity OFF
doc = SettingsDocument.from_file(src); doc.validate(), doc.get("useGravity")   # [] , True
doc.save(out); json.loads(out.read_text())["useGravity"]               # true
json_to_config({"useGravity": True}).useGravity is True                # True   -> reducer: gravity ON
```

At the base 7b6d6b9 the same file reports `Gravity correction (useGravity): expected true/false, got 1` and saves `1`
(reducer reading unchanged). The design reviewer reproduced the same, independently, plus the injected path
(`SettingsDocument(config)` with `useGravity=1` → `validate() == []`).

**Why it blocks.** Load → save in the editor silently changes a reduction parameter, with no problem line, where the
base flagged it. Reachability: the plan's own writer list includes "a hand-edited file (anything)", and the scientists'
item 1 says they spell booleans as 1/0. Frequency, stated honestly: **0 of 79** real `*_settings.json` /
`reduce_settings.json` under `IPTS-3[6-8]*` hold an int here (79 `true`, 6 absent) — reachable, not yet observed.

**The premise that failed.** Plan §9 / learning §4 "0/1 is what the consumer writes and reads". True for `useBS`
(`:509`, `:979`) and for six of the seven scalar booleans; false for `useGravity`:

| Field | Reader | Reads by |
|---|---|---|
| `useBS[i]` | `nr_reduction_calc.py:509`, `:979` | truthiness |
| `Normalize` | `:1137` | truthiness |
| `AutoScale` | `:156` | truthiness |
| `plotON` | `:135`, `:236`, `:474`, `:584` | truthiness |
| `plotQ4` | `:223` | truthiness |
| `save8col` | `new_reduction_from_file.py:261` (and `nr_reduction_calc.py:128-129`) | truthiness |
| `use_emission_time` | `:430`, `:438` | truthiness |
| **`useGravity`** | **`:1079`** | **identity (`is True`)** |

## Fix — behaviour, not code (amendment 18)

**Type domain.** 8 declared booleans: 7 scalars + 1 list (`useBS`). Mixed consumers (7 truthiness, 1 identity), so the
prescription is behavioural. The reproduction covered **one scalar (`useGravity`) with int `1`**. Int `0` is harmless
there (held `False`, saved `false`, reducer off both times) but must still be in the guard.

**Required behaviour (B8, new):** *load → save never changes what the reducer does with a declared boolean.* For a
field the reducer reads by truthiness, 1/0 is accepted, canonicalized and saved exactly as v1 does. For a field read by
identity (`useGravity` today), an int keeps the base behaviour: reported (the message may name the accepted spelling),
held as written, saved as written. How the distinction is declared (a `Field` attribute; a narrower `as_boolean` for
that field; …) is the plan owner's choice. One declaration, no name test in the document; the v1 single-dispatch design
otherwise stands.

**Guard (must vary the dimension the fix freezes).** Parametrize over **every declared scalar boolean × {1, 0}** and
`useBS` entries × {1, 0}. Assert: the value the reducer acts on (its reading per the table above) is the same for the
source file and for the editor-saved file. Assert on `type()`/saved text, not `==`. A guard written for `useGravity`
alone freezes the field dimension and would miss the next identity reader. Mutation row: drop the identity-field
distinction → the guard reds on `useGravity`=1 only.

**Out of scope here, for the Analyst/human:** changing `:1079` to truthiness would let v1's A3 stand. It is a reducer
file (OUT) and a numerical-behaviour change (`# TODO: Implementation needs checking/deciding whether to keep!` on that
line), so it belongs to a separate slug or the scientists' ruling, not this one.

## Advisories (non-blocking; carried to the PR body)

Design reviewer:
- `settings_editor.py:248` scalar checkbox renders via `bool(value)` (pre-existing), a second boolean interpretation
  outside `as_boolean`: `"0"`/`2` show ticked. The "one place that decides" docstring and commit claim overstate this.
- `settings_document.py:344`, `:392` repeat `self._encode_for_file(make_json_safe(self.to_dict()))`. A `_file_payload()`
  would make "one encoder" structural.
- `settings_document.py:292` / `settings_editor.py:206`: the record exemption is keyed on `runtime_owned` (not
  per-angle). Today that means exactly the two Lambda fields; a future scalar runtime-owned field would be silently
  exempt. Consider `Field.record_only` or a comment.
- File sizes over the 300-line guideline (`field_spec.py` 632, `settings_editor.py` 481, `settings_document.py` 437;
  pre-existing). Queue a split.

Test reviewer:
- V4 asserts the `LambdaMinUse` display only, not `LambdaMaxUse`.
- V5 checks that `useBS` appears in the report, not that it appears once (§5 row).
- No view test for §5 "scalar `Normalize: 1` → checkbox checked" or for loading a non-list `useBS` through the tab.
- `_build_editor`'s runtime-owned branch calls `_show` redundantly (`refresh_scalars` re-shows); removing it changes
  nothing (176 passed under the mutant).

Integrator:
- §8.5 names `shared/autoreduce/` folders. Reducer-written `*_settings.json` live under `shared/reduced/`.
- `Aug2026/REFL_231105_settings.json` fails the acceptance script only on a per-angle **count** line naming `useBS`,
  identical at base and tip (the script filters by field name) → ledger `todo-real-settings-per-angle-length-mismatch.md`.

— Integrator, Claude Opus 5.5
