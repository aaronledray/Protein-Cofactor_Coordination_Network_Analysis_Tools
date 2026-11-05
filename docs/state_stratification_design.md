# State-stratified conserved-contact analysis (design)

## Purpose

Show how residue contacts vary across explicitly labelled structural states in a sequence set, while retaining the existing pooled sequence vote as the primary result. This is a descriptive comparison of the structures supplied; it does not establish a causal state effect or a population frequency.

## Interface

- Add an opt-in `sscna conserve --stratify-state` flag.
- Require a non-empty `state` value for every manifest row. State labels are global and case-sensitive: the same label must mean the same state in every sequence and group.
- Leave the current pooled `conserved_contacts` and all existing defaults unchanged.
- Add a tidy `state_contacts` table to the same log, with one row per contact candidate and state. Preserve the state label, per-group positive and total sequence votes, support fractions, and whether each group passes the existing support and minimum-sequence rules.
- Add a `state_structure_votes` table with one row per candidate pair, state, group, and sequence. Record eligible structures, positive structures, their fraction, and the resulting sequence vote; a sequence-state with no eligible structures has a missing vote rather than a negative.
- Add a `state_contrasts` table with one row per candidate pair and unordered pair of states. Report each state's overall status (`pass`, `support_fail`, or `insufficient_coverage`) and a comparison category. This is a descriptive summary of the existing vote results.
- Keep audit behavior unchanged. State-specific directives are out of scope for the first version; reject `--directives` with `--stratify-state` because existing directives describe only the pooled vote.

## Vote and denominators

1. Build the candidate contact set as the union of residue pairs observed in contact in any supplied state, using the existing residue-contact rules and cutoff.
2. Within each `(pair, sequence_id, state)`, count only structures in that state where both positions are mapped. A structure missing either position is not a negative. The sequence-state vote is positive when at least half of its eligible structures show the contact, matching the existing `structure_fraction=0.5` rule.
3. For each `(pair, state, group)`, count sequence-state votes. Apply the existing defaults: at least five voting sequences and at least 90% positive in every group. A sequence with no eligible structure for that pair and state contributes no vote.
4. Report each state's voter count and support fraction, including states that fail the minimum or support threshold. Do not collapse a no-voter state into a negative result; mark it as insufficient coverage. Preserve per-sequence eligible and positive structure counts in `state_structure_votes` so repeated structures remain auditable.
5. Vote chemistry flags only among contact-positive structures within each sequence-state, matching the existing conditional chemistry vote. The existing opt-in sulfur-acceptor and sequence-neighbor rules also apply when requested.

## Interpretation and limits

- State-stratified votes answer whether a contact meets the existing conservation rule within each labelled state. They do not test whether states differ statistically.
- Structures remain equally weighted inside each sequence-state vote. Multiple structures from the same crystal form or experimental setup can therefore dominate that within-state vote; report raw structure counts and retain the current non-independence caveat.
- State labels and residue correspondence remain caller-supplied. `anchor_audit` and `backbone_only_audit` remain separate input-QC tools; stratification does not silently drop structures or sequences.
- Cofactor, water-bridge, and relay state stratification is not included in the first version. The initial output covers residue-residue contacts only.
- Pairwise comparison labels summarize vote status only; they do not estimate a state effect or its uncertainty.

## Validation

- Confirm that without `--stratify-state`, the log and pooled vote are unchanged, including with a manifest that has a `state` column.
- Use a small hand-built fixture with two states to check positive, negative, missing-position, no-voter, and below-minimum cases, including a sequence represented in only one state.
- Check that the per-state vote equals running the existing vote independently on each state's subset, with the union candidate set retained so state-specific absences are visible.
- Check that state-specific support and chemistry are not pooled across states, and that insufficient coverage is explicit in both CSV and summary.
- Keep the existing oxidase 490-contact reproduction and stable-core schemas unchanged.

## Contrast categories

The pairwise summary distinguishes `pass_both`, `pass_only_a`, `pass_only_b`, support failure in both states, insufficient coverage in both states, and the combinations where one state passes or fails while the other has insufficient coverage. State A and B are the sorted labels in the output row. A pass requires every group to meet its support and minimum-sequence rules; insufficient coverage takes precedence over support failure when any group is below the minimum.
