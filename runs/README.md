# runs

Output from `examples/GbFacetMeshing`, written here rather than into a build directory so that
rebuilding, reconfiguring or clearing a build tree cannot touch results.

Layout, one directory per boundary and one per set of sweep settings inside it:

    runs/sigma<N>_theta<deg>_axis<h.k.l>_plane<h.k.l>/<search>_maxEngaged<N>/

Each run directory holds:

| file | covers |
|---|---|
| `states.txt` | the states that were built: energies, coincidence counts, why each was chosen |
| `signatures.txt` | their full signatures, in the layout `fin` reads back |
| `output_thread_*.txt` | one line per state **surveyed**, with its full signature |
| `sortedByUnrelaxed/Tethered/Full.txt` | every surveyed state, ordered by each energy |
| `state_<i>_0.txt` / `_1.txt` | undeformed and deformed configurations, extended XYZ |
| `dump.state_<i>_2` / `_3` | tethered and freely relaxed structures, LAMMPS dumps |
| `shortlist/engaged_<n>.xyz` | one trajectory per engagement level, for viewing |

`shortlistStates.py` builds the `shortlist/` trajectories; run it with `--run-dir` pointing at a
run directory, and `--top` matching the sweep's `shortlistTop`.
