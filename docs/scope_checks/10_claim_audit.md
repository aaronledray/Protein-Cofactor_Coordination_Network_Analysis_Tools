# Claim audit (scope check 5)

The question: does anything in the repository's text, figures or generated reports say more than the evidence supports? Method: read the tracked docs, module and
script docstrings, the legends that `class_rules` generates, figure captions and the generated scope-check reports for strong or predictive wording
(important, critical, stabilize, design, ensures, textbook, real, validated, exactly, no hydrogens), and compare each claim with what was measured.

## Found and fixed

1. **No statement of what the output means.** The docs said what the module computes but not what a conserved contact does and does not imply. Added "What the output
   means" to `docs/conserved_interactions.md`: the output is a description of recurrence and chemistry; it does not say a contact is important or that changing it matters;
   conservation can reflect ancestry; the project's TEM-1 result (outside this repository and not rerun here) was negative for predictive use of contact-derived directives, so a
   directive should be worded as a description with group sizes and margins beside it.
2. **"No hydrogens are used" was out of date.** The opt-in flags `polar_s_geom` and the neighbor `polar` class place a hydrogen from the backbone or ring. The module
   docstring, the `polar` legend text and the docs now say hydrogens are placed only for those flags and that a placed hydrogen is a model, not an observation.
3. **A circular reading of the plastocyanin result.** "The vote recovers the copper site" was true but not independent: the Cu tier cutoffs come from a survey that includes
   plastocyanins. The plastocyanin write-up and the generated vote report now say the recovery is a consistency check, not an independent test.
4. **An unsourced functional claim.** An earlier draft attached the acidic patch of plastocyanin to a retracted contact. It had no source in this repository; the sentence
   was removed and the retraction stated plainly.
5. **A figure caption asserted a hydrogen bond.** Figure 3 said a backbone N-H within 3.7 A "is a hydrogen-bond donor". It now says "candidate" and points to the placed-hydrogen
   geometry shown on the figure (the survey: about nine in ten such contacts to a metal-bound Cys point the N-H at the sulfur, about half to a free Cys).
6. **Overreach in the shuffle report.** "The vote does not manufacture conservation" held for one null only. The report now says the null is a relabelling of positions
   within each structure, that other nulls (unrelated proteins, a different alignment) were not tested, and that the oxidase rows use three replicates and depend on which groups the scrambled
   sequences fall in.
7. **Neighbor polar wording.** The neighbor `polar` class was called "a hydrogen bond"; for hydroxyl and Lys donors the hydrogen is not placed and the distance decides, so it is now
   "a hydrogen-bond candidate".

## Checked and left alone

- The 490-contact reproduction is described as using the manuscript's own residue map and as not testing `load_atoms` (already stated in the Validation section).
- `polar`, `polar_s*` and the neighbor classes are described as capable pairs or candidates throughout; no text calls a flag a proven interaction.
- The Met tiers (2.7 and 3.3 A) are labelled conventions wherever they appear.
- The surveys state their samples (atomic-resolution X-ray, one chain per 30% cluster, first model) and that hydrogens are placed.
- The rules legend (`class_rules`) is generated from the parameters used and is checked against the code at every documented cutoff.
- Directives (added 2026-11-02) are restatements of table rows with their evidence, generated from fixed templates; the tests refuse a directive with another tier or claim level and scan the statements, the policy and everything the log writes for importance or predictive wording.

## Not resolved here

- The TEM-1 statement rests on the project notes, not on anything in this repository; it should be cited to the manuscript or design log when the log format is chosen.
- Whether a design log may use the word "conserved" without the group sizes and margins is the author's call (the stability section recommends carrying them).
