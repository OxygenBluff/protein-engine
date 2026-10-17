# AAPFA: Peptide Sequence Generator

Generates peptide sequences that satisfy user-defined physicochemical
constraints, and lets you watch the search algorithm work, step by step.
Built with classical AI (no machine learning), as a homage to my Biology
baccalaureate roots.


## Inputs

```python
targets = {
    'hydro_min': 8.25,   'hydro_max': 8.35,   # hydrophobicity range
    'mass_min': 10500,   'mass_max': 11500,   # molecular mass range
    'stability_min': 80.0,   # stability = based on the aliphatic index
    'binding_min': 220.0,    # binding affinity = polarity + |charge|
}
target_length = 100
biological_switch = True
```

**Biological switch.** When on, candidates must also be biologically
plausible: a sliding-window folding check and charged N/C-terminus (hooks)
checks.
When off, the problem becomes purely mathematical: any sequence
meeting the numeric constraints is valid, whether or not it could exist
in nature.

## How it works

Sequences are built by a **depth-first search (DFS) with aggressive pruning**,
guided by a **greedy heuristic**.

**Greedy ordering.** At each step, candidate amino acids are ranked by how
well they move the sequence toward the targets, with a penalty for repeating
the same residue. The heuristic decides which branch DFS tries first. It never
decides what is valid. 

**Pruning.** The search space explodes with length, so a branch is cut the
moment it can no longer satisfy a constraint.

**Branch and bound ** at each position, the algorithm computes
the best and worst case value of every property over the remaining slots
(hydrophobicity, mass, stability, binding affinity). If the target range is
unreachable even in the best case, the whole subtree is discarded unexplored.
For hydrophobicity, stability and binding the bounds are exact (per-residue
extremes). For mass, I use a realistic residue range (100-150 Da) instead
of the theoretical extremes, which would be too loose.

**Explainable by design:** every pruning decision the search makes is logged with its exact reason (which constraint, which bound,which residue), so any result can be traced back step by step.


```
109 | PRUNE | 25 | Phe | molecular mass exception: Too heavy | Val-Pro-Arg-Ser-...-Phe
```

Columns: step | event | depth | amino acid | reason | current sequence.
The full log can be downloaded once generation finishes.

## Features

- Sequence generation from numeric constraints, with optional biological realism
- Step-by-step visualization of the search, including every prune and why
- Downloadable log file
- History table of generated sequences 
- **Compare any 2 sequences** (generated or typed in by hand) using
  per-residue **profiles** (feature value vs. position) and similarity
  scores 

## Architecture

| Layer | Tech | Role |
|-------|------|------|
| Core algorithm | Python (`AA_Math`, `AA_Validation`, `AA_BD`, `AA_API`) | Generation, validation, scoring |
| Backend | Spring Boot (Java), Streams the very large step-by-step log |
| Frontend | Angular | Visualization and comparison screens |

The log is massive, so everything is asynchronous and streamed instead of
returned in one response. 

## Running it

**Requirements:** Python 3.x. The generation core (`AA_Math.py`) uses only the
standard library. The API and comparison features need:
	pip install -r requirements.txt 

**Python core:** edit the targets (`targets`, `target_length`,
`biological_switch`) at the bottom of `AA_Math.py`, then run:
	python AA_math.py 

**Backend**
```bash
cd java_backend
[./mvnw spring-boot:run]
```

**Frontend**
```bash
cd frontend
npm install
ng serve
```
Then open `http://localhost:4200/`.

## Limitations & future work

- Classical AI only: no learned model.
- Valid, not optimal: the first sequence that passes all constraints is
  returned.
- The Java backend works but could use a refactor.
- Idea: an ML model to score or guide candidates, replacing some of the
  hand-written pruning rules.

## Author

Mhamdi Mohamed Habib,https://www.linkedin.com/in/mhamdi-mohamed-habib-4787433b3