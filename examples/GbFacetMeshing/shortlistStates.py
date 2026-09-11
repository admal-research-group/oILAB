#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Collect the mesostates of a GbFacetMeshing sweep into one trajectory per level.

The sweep itself now does the shortlisting: it surveys every state, keeps only the
numbers, and rebuilds the shortlisted states alone.  What lands in its output directory
is therefore already the shortlist -- states.txt plus, per state, four files:

    state_<index>_0.txt   undeformed,           extended XYZ
    state_<index>_1.txt   deformed,             extended XYZ
    dump.state_<index>_2  tethered relaxation,  LAMMPS dump
    dump.state_<index>_3  free relaxation,      LAMMPS dump

beside signatures.txt, which holds the full signature of each shortlisted state in the
layout the sweep reads signatures back in -- point `fin` at it with enumerateStates
false and those states are rebuilt.  Every state's signature, shortlisted or not, is in
the sweep's output_thread_<id>.txt files.  Nothing here reads either; they are what
makes a state in these trajectories traceable back to the mesostate it came from.

This script's job is to turn that into something a viewer can open.  It re-derives the
same ranking, which on a directory the sweep already shortlisted is a check rather than
a filter -- run it with a smaller --top to narrow further, or on an older single-pass
run to shortlist one for the first time.

States that engage different numbers of coincidence points are not comparable -- a
state holding eight nodes costs more than one holding a single node whatever the
boundary does -- so the ranking is taken within each level of engagement separately:
the lowest `--top` states at one coincidence, the lowest `--top` at two, and so on to
the ceiling the sweep reached.  The level is `nodes`, the number of coincidences the
deformed structure holds -- not the number the enumeration engaged, which is `engaged`
and describes the route rather than the boundary.

Within a level the enumeration produces whole families of states that differ only by
where along the periodic directions their engaged nodes sit -- translations and
symmetry images of one configuration -- and these are degenerate: every number the
sweep records for them, energy and density and the relaxation figures alike, agrees to
the last digit written.  Left alone one such family fills the whole top `--top` with
copies of a single structure and hides the next distinct one, so states agreeing on
all of those numbers are treated as one and only the first is kept, the rest counted
against it in the `copies` column of the summary.  `--no-dedup` restores the plain
ranking.

Each level's shortlist is written as a single extended-XYZ trajectory,
`engaged_<nodes>.xyz`, holding its kept states end to end in rank order, four frames
per state: undeformed, deformed, tethered relaxation, free relaxation, in that order.  A
viewer opens one such file and has the whole level as one sequence -- state 1 at frames
0-3, state 2 at frames 4-7, and so on -- so both the construction -> relaxation of a
single state and the progression from best to worst can be stepped through without
loading a file per state and keeping them in order by hand.  `--per-state` writes the
individual per-state files as well, for pulling one state out on its own.

The two relaxations are separate LAMMPS runs from the same starting configuration: the
tethered one holds the boundary atoms where the construction put them, so it costs the
state the enumeration built; the free one lets them go, so it costs the boundary that
state turns into.  A run with no tether writes only the first, and those states come
out as three frames.

Every frame repeats its state's identity and energies as extended-XYZ comment keys
(`state`, `rank`, `nodes`, `energy`, ... and `stage`/`Frame` for which of the three
configurations it is), which OVITO exposes as global attributes: a Text label modifier
reading `[state]` / `[Frame]` / `[energy]` then names on screen whatever frame is
showing, which is what makes a sixty-frame file navigable.

The atom count changes between frames -- LAMMPS deletes atoms that the deformation
brings on top of one another -- which extended XYZ carries without complaint, each
frame declaring its own count.  How many it deletes is the `fused` column, and it is
generally more than the state engaged: the displacement field cannot be aimed at the
engaged nodes alone, so it closes other pairs of atoms along with them.

Usage:
    ./shortlistStates.py                       # ./generalGB_* found automatically
    ./shortlistStates.py --run-dir <dir> --top 20
    ./shortlistStates.py --sort-by full
"""

import argparse
import glob
import os
import sys

# The columns states.txt may carry, in the order the sweep writes them.  Which of the
# optional ones are present depends on how the run was configured -- unrelaxed runs
# have no `beforeRelaxation`, untethered ones no `spring` -- so the header line is
# read rather than assumed.
REQUIRED_COLUMNS = ("index", "nodes")

# Columns that say how a state was reached, or how many others were folded into it, rather than
# what it is; they are left out of the fingerprint that decides whether two states are the same.
# `engaged` is the number of nodes the enumeration engaged; `nodes` is the number of coincidences
# the structure ended up with, and the sweep can reach one structure by engaging different numbers
# of its coincidences -- the rest forming on their own -- so the same boundary arrives several
# times over with `engaged` the only column separating the copies.  Keying on it would preserve
# exactly the duplicates the dedup exists to remove.  `copies` is the sweep's own count of the
# family it already folded, which says nothing about the structure.
PATH_COLUMNS = ("engaged", "copies")


def stateFingerprint(state, header):
    """The numbers that have to agree for two states to be one and the same.

    Distinct placements of the engaged nodes that are translations or symmetry images
    of one another relax to the same structure, and the sweep has no column that tells
    them apart -- the energies, the density and the corrugation all agree.  Every
    numeric column is included rather than the ranked one alone: two states really are
    the same only if nothing recorded separates them, and a bare energy match would
    merge an accidental degeneracy between genuinely different structures.

    PATH_COLUMNS are excluded: they record the route the enumeration took rather than
    the structure it arrived at, and two states differing only there are one boundary.
    """
    return tuple(state[column] for column in header
                 if column not in ("index", "nodes")
                 and column not in PATH_COLUMNS
                 and isinstance(state.get(column), float))


def groupByFingerprint(states, header, tolerance):
    """A group index per state, equal for states whose fingerprints agree within `tolerance`.

    This mirrors groupByFingerprint() in GbFacetMeshing.cpp, so that re-deriving the
    sweep's own shortlist here reproduces it rather than disagreeing with it at the
    margin.  The comparison is a tolerance rather than a rounded key: a rounded key is
    cheaper but splits two values falling either side of a bucket edge, and that is
    exactly what a minimiser's last-digit noise does -- two copies of one structure
    whose free relaxations came back 4e-8 apart land in different buckets at any
    tolerance whose edge lies between them.  Sorting on the fingerprint puts such
    states next to one another, and each is compared against its group's first member
    rather than its predecessor so that a run of near-equal values cannot drift a
    group away from where it started.
    """
    fingerprints = [stateFingerprint(state, header) for state in states]
    order = sorted(range(len(states)), key=lambda i: fingerprints[i])
    group = [None] * len(states)
    representative, groups = None, 0
    for i in order:
        matches = representative is not None and all(
            abs(mine - theirs) <= tolerance
            for mine, theirs in zip(fingerprints[i], fingerprints[representative]))
        if not matches:
            representative, groups = i, groups + 1
        group[i] = groups - 1
    return group


def parseManifest(path):
    """Read states.txt into a list of dicts, one per state.

    The trailing `engaged nodes` field is free text describing each engaged node, so
    only the leading numeric columns named in the header are parsed; the description
    is kept whole under 'description'.
    """
    with open(path) as manifest:
        lines = manifest.readlines()

    header = None
    for line in lines:
        if line.startswith("#") and "index" in line and "nodes" in line:
            header = line.lstrip("#").split()
            break
    if header is None:
        raise RuntimeError("no column header found in " + path)

    # `engaged nodes` is two words naming one trailing free-text field.
    if header[-2:] == ["engaged", "nodes"]:
        header = header[:-2]
    for column in REQUIRED_COLUMNS:
        if column not in header:
            raise RuntimeError("column '%s' missing from the header of %s" % (column, path))

    states = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        if len(fields) < len(header):
            continue
        state = {"description": " ".join(fields[len(header):])}
        for column, value in zip(header, fields[:len(header)]):
            state[column] = value
        state["nodes"] = int(state["nodes"])
        for column in header:
            if column not in ("index", "nodes"):
                try:
                    state[column] = float(state[column])
                except ValueError:
                    pass
        states.append(state)
    return states, header


def readExtendedXyz(path):
    """(count, commentLine, [atom lines]) from an extended XYZ file."""
    with open(path) as configuration:
        count = int(configuration.readline().split()[0])
        comment = configuration.readline().rstrip("\n")
        atoms = [configuration.readline().rstrip("\n") for _ in range(count)]
    return count, comment, atoms


def latticeFromDumpBounds(bounds, tilts):
    """The cell vectors and origin of a LAMMPS dump box.

    A dump records the *bounding box* of a triclinic cell, which is the cell grown to
    contain its own tilt; the tilts have to be taken back out to recover the edges.
    For an orthogonal box the tilts are zero and this is the identity.

    The origin matters as much as the edges here: LAMMPS pads the non-periodic
    direction, so its box neither starts nor ends where the extended-XYZ frames' does,
    and a frame given edges without an origin draws its cell away from its atoms.
    """
    (xloBound, xhiBound), (yloBound, yhiBound), (zloBound, zhiBound) = bounds
    xy, xz, yz = tilts
    xlo = xloBound - min(0.0, xy, xz, xy + xz)
    xhi = xhiBound - max(0.0, xy, xz, xy + xz)
    ylo = yloBound - min(0.0, yz)
    yhi = yhiBound - max(0.0, yz)
    lattice = [[xhi - xlo, 0.0, 0.0],
               [xy, yhi - ylo, 0.0],
               [xz, yz, zhiBound - zloBound]]
    return lattice, [xlo, ylo, zloBound]


def readLammpsDump(path):
    """(count, lattice, origin, [(id, type, x, y, z)]) from a LAMMPS dump's last frame."""
    with open(path) as dump:
        lines = dump.readlines()

    # The minimization writes one frame, but a restarted run can leave several; the
    # last is the converged one.
    starts = [i for i, line in enumerate(lines) if line.startswith("ITEM: TIMESTEP")]
    if not starts:
        raise RuntimeError("no frame found in " + path)
    frame = lines[starts[-1]:]

    count = None
    bounds, tilts = [], [0.0, 0.0, 0.0]
    atoms = []
    index = 0
    while index < len(frame):
        item = frame[index]
        if item.startswith("ITEM: NUMBER OF ATOMS"):
            count = int(frame[index + 1])
            index += 2
        elif item.startswith("ITEM: BOX BOUNDS"):
            triclinic = "xy" in item
            bounds, tilts = [], [0.0, 0.0, 0.0]
            for axis in range(3):
                values = [float(value) for value in frame[index + 1 + axis].split()]
                bounds.append((values[0], values[1]))
                if triclinic and len(values) > 2:
                    tilts[axis] = values[2]
            index += 4
        elif item.startswith("ITEM: ATOMS"):
            columns = item.split()[2:]
            required = ("id", "type", "x", "y", "z")
            if not all(column in columns for column in required):
                raise RuntimeError("dump %s lacks one of %s" % (path, str(required)))
            where = {column: columns.index(column) for column in required}
            for line in frame[index + 1: index + 1 + (count or 0)]:
                fields = line.split()
                if len(fields) < len(columns):
                    break
                atoms.append((int(fields[where["id"]]),
                              int(fields[where["type"]]),
                              float(fields[where["x"]]),
                              float(fields[where["y"]]),
                              float(fields[where["z"]])))
            index += 1 + len(atoms)
        else:
            index += 1

    # LAMMPS emits atoms in whatever order its domain decomposition holds them.
    atoms.sort(key=lambda atom: atom[0])
    lattice, origin = latticeFromDumpBounds(bounds, tilts)
    return len(atoms), lattice, origin, atoms


def commentFields(comment):
    """The key="value" and key=value pairs of an extended XYZ comment line."""
    fields, key, value, quoted, reading = {}, "", "", False, False
    token = ""
    for character in comment:
        if character == '"':
            quoted = not quoted
            continue
        if character == "=" and not quoted and not reading:
            key, token, reading = token.strip(), "", True
            continue
        if character == " " and not quoted and reading:
            fields[key] = token.strip()
            key, token, reading = "", "", False
            continue
        token += character
    if reading:
        fields[key] = token.strip()
    return fields


def appendFrames(trajectory, frames, state, rank, sortColumn, firstTime):
    """Append one state's configurations to an open extended-XYZ trajectory.

    Each frame keeps the Lattice and Properties of the file it came from, so the
    frames stay readable exactly as the sweep's own files are, and gains the state's
    identity and energy.

    `Time` counts from `firstTime` rather than from zero so that many states can share
    one file and still have each frame's Time match the frame number the viewer shows.
    Which of the three configurations a frame holds is then carried by `stage` and
    `Frame`, not by its position: in a level file frame 37 is the deformed
    configuration of the thirteenth state, and only those keys say so.

    Returns the Time the next state should start at.
    """
    for offset, frame in enumerate(frames):
        label, count, lattice, properties, atoms, extra = frame
        trajectory.write("%d\n" % count)
        comment = ['Lattice="%s"' % lattice,
                   "Properties=%s" % properties]
        # PBC and origin come from the frame's own source: LAMMPS pads the
        # non-periodic direction, so the minimized frame's box is not the
        # constructed one and must not be given its origin.
        for key in ("PBC", "origin"):
            if extra.get(key):
                comment.append('%s="%s"' % (key, extra[key]))
        comment += ["Time=%d" % (firstTime + offset),
                    "stage=%d" % offset,
                    'Frame="%s"' % label,
                    'state="%s"' % state["index"],
                    "rank=%d" % rank,
                    "nodes=%d" % state["nodes"]]
        for column in ("density", "unrelaxed", "tethered", "spring", "full",
                       "corrugation"):
            if column in state and isinstance(state[column], float):
                comment.append("%s=%.8f" % (column, state[column]))
        # These are counts; eight decimal places on them would be noise dressed as detail.
        for column in ("engaged", "fused", "expelled", "copies"):
            if column in state and isinstance(state[column], float):
                comment.append("%s=%d" % (column, int(round(state[column]))))
        comment.append('rankedBy="%s"' % sortColumn)
        trajectory.write(" ".join(comment) + "\n")
        for atom in atoms:
            trajectory.write(atom + "\n")
    return firstTime + len(frames)


def collectState(runDirectory, state):
    """Read one state's three configurations.

    Returns the list of frames that were found; a state whose minimized dump is
    missing -- an unrelaxed run, or a LAMMPS failure -- still gets its first two.
    Nothing is written here: the same frames go to the level's trajectory and,
    optionally, to the state's own file, and reading the dump twice for that would
    double the cost of the whole shortlist.
    """
    index = state["index"]
    frames = []

    for suffix, label in ((0, "undeformed"), (1, "deformed")):
        path = os.path.join(runDirectory, "state_%s_%d.txt" % (index, suffix))
        if not os.path.exists(path):
            return frames, "missing " + os.path.basename(path)
        count, comment, atoms = readExtendedXyz(path)
        fields = commentFields(comment)
        frames.append((label, count,
                       fields.get("Lattice", ""),
                       fields.get("Properties", "atom_types:I:1:pos:R:3:radius:R:1"),
                       atoms,
                       {"PBC": fields.get("PBC", "F T T"),
                        "origin": fields.get("origin", "")}))

    # The dump carries no radius, so the deformed frame's is reused: it is a drawing
    # hint, and giving a relaxed frame a different one would make the same atom change
    # size between frames for no physical reason.
    radius = "0.05"
    if frames[-1][4]:
        deformedFields = frames[-1][4][0].split()
        if len(deformedFields) >= 5:
            radius = deformedFields[4]
    deformedPbc = frames[-1][5].get("PBC", "F T T")

    # The sweep relaxes each state twice from the same starting configuration and keeps
    # both structures: _2 held the boundary atoms where the construction put them, _3
    # let them go.  They are different questions about the same state -- what the
    # constructed boundary costs, and what boundary it turns into -- so both are frames
    # of the trajectory, in that order, and a run with no tether writes only _2.
    for suffix, label in ((2, "tethered"), (3, "relaxed")):
        dumpPath = os.path.join(runDirectory, "dump.state_%s_%d" % (index, suffix))
        if not os.path.exists(dumpPath):
            continue
        count, lattice, origin, atoms = readLammpsDump(dumpPath)
        lines = ["%d %.8f %.8f %.8f %s" % (atom[1], atom[2], atom[3], atom[4], radius)
                 for atom in atoms]
        frames.append((label, count,
                       " ".join("%.15g" % value for row in lattice for value in row),
                       "atom_types:I:1:pos:R:3:radius:R:1", lines,
                       {"PBC": deformedPbc,
                        "origin": " ".join("%.15g" % value for value in origin)}))

    return frames, None


def main():
    parser = argparse.ArgumentParser(
        description="Shortlist the lowest-energy states of a GbFacetMeshing sweep, "
                    "one ranking per level of engagement, and write each level as a "
                    "single trajectory holding its states' undeformed, deformed and "
                    "minimized configurations as consecutive frames.")
    parser.add_argument("--run-dir", default=None,
                        help="the sweep's output directory (default: the only "
                             "generalGB_* directory here)")
    parser.add_argument("--out-dir", default=None,
                        help="where the shortlist is written "
                             "(default: <run-dir>/shortlist)")
    parser.add_argument("--top", type=int, default=20,
                        help="how many states to keep per level (default 20)")
    parser.add_argument("--sort-by", default="tethered",
                        help="the column ranked on: tethered (default; the relaxation that "
                             "holds the boundary atoms where the construction put them), full "
                             "(the free relaxation), or unrelaxed")
    parser.add_argument("--no-dedup", action="store_true",
                        help="keep every state, including ones the sweep records "
                             "identically (translations and symmetry images of the "
                             "same structure)")
    parser.add_argument("--redo-dedup", action="store_true",
                        help="fold duplicates again even on a directory the sweep has "
                             "already deduplicated; see the note in the source on why "
                             "that can disagree with the sweep")
    parser.add_argument("--dedup-tol", type=float, default=1e-6,
                        help="how close two numeric columns count as equal when "
                             "deduplicating (default 1e-6)")
    parser.add_argument("--per-state", action="store_true",
                        help="also write each state to its own file, beside the "
                             "level trajectory that holds them all")
    parser.add_argument("--flat", action="store_true",
                        help="with --per-state, write every state file into one "
                             "directory instead of one subdirectory per level")
    arguments = parser.parse_args()

    runDirectory = arguments.run_dir
    if runDirectory is None:
        # Any directory holding a states.txt is a run, whatever it is called.  The sweep names
        # its output for the boundary -- sigma, misorientation, axis, plane -- with one folder
        # per set of run settings inside that, so there is no fixed prefix to match on.
        candidates = sorted(os.path.dirname(path)
                            for path in glob.glob("*/states.txt") + glob.glob("*/*/states.txt"))
        if len(candidates) != 1:
            parser.error("found %d run directory(s) below here%s; pass --run-dir"
                         % (len(candidates),
                            "" if not candidates else ":\n  " + "\n  ".join(candidates)))
        runDirectory = candidates[0]
    manifestPath = os.path.join(runDirectory, "states.txt")
    if not os.path.exists(manifestPath):
        parser.error("no states.txt in " + runDirectory)

    outputDirectory = arguments.out_dir or os.path.join(runDirectory, "shortlist")
    os.makedirs(outputDirectory, exist_ok=True)

    states, header = parseManifest(manifestPath)
    # Whether to fold duplicates here at all.
    #
    # A `copies` column means the sweep already did it, and did it better: it grouped each
    # state against the whole level -- hundreds of thousands of states -- while all this
    # script has is the few dozen the sweep kept.  Grouping by tolerance is not stable
    # under taking a subset, because which state becomes a group's representative depends
    # on which states are present: two states the sweep separated by comparing both
    # against a third can fall within tolerance of each other once that third is gone.
    # Re-folding here would therefore quietly contradict the sweep's own count rather
    # than confirm it, so on a two-pass directory the sweep's answer stands and `copies`
    # is read from the manifest.  --redo-dedup overrides this; --no-dedup on an older
    # single-pass run turns folding off altogether.
    sweepDeduplicated = "copies" in header
    if not states:
        parser.error("states.txt holds no states")
    if arguments.sort_by not in header:
        parser.error("no column '%s' in states.txt; available: %s"
                     % (arguments.sort_by, ", ".join(header)))
    if not isinstance(states[0].get(arguments.sort_by), float):
        parser.error("column '%s' is not numeric -- was the sweep run with energies "
                     "off?" % arguments.sort_by)

    print("read %d states from %s" % (len(states), manifestPath))
    dedup = not arguments.no_dedup and (arguments.redo_dedup or not sweepDeduplicated)
    if sweepDeduplicated and not arguments.redo_dedup and not arguments.no_dedup:
        print("the sweep already folded duplicates; its `copies` column is used as it stands")

    byLevel = {}
    for state in states:
        byLevel.setdefault(state["nodes"], []).append(state)

    summaryPath = os.path.join(outputDirectory, "shortlist.txt")
    summary = open(summaryPath, "w")
    summary.write("# lowest %d states per level of engagement, ranked by %s\n"
                  % (arguments.top, arguments.sort_by))
    # `copies` is appended in a column of its own -- read from the manifest when the sweep
    # deduplicated, counted here when it did not -- so it is dropped from the passed-through
    # columns rather than printed twice.
    passedThrough = [column for column in header
                     if column not in ("index", "nodes", "copies")]
    summary.write("# nodes  rank  index  %s  copies  frames\n"
                  % "  ".join(passedThrough))
    summary.write("# `frames` locates the state in its level trajectory as "
                  "file:first-last: undeformed, deformed, then the tethered and the "
                  "freely relaxed structures\n")

    written, incomplete, levelFiles = 0, [], []
    for level in sorted(byLevel):
        ranked = sorted(byLevel[level], key=lambda state: state[arguments.sort_by])

        # The whole level is scanned even once `--top` distinct states are in hand:
        # the copies of a kept state are spread through the ranking and counting them
        # is what turns "20 of 258" into a statement about how large each family is.
        # When the sweep deduplicated, its own count of each family is what gets reported.
        copies = {state["index"]: int(state["copies"]) for state in ranked} \
            if sweepDeduplicated else {}
        if not dedup:
            keep = ranked[:arguments.top]
        else:
            group = groupByFingerprint(ranked, header, arguments.dedup_tol)
            keep, representative = [], {}
            for position, state in enumerate(ranked):
                key = group[position]
                if key in representative:
                    copies[representative[key]] += 1
                    continue
                if len(keep) < arguments.top:
                    representative[key] = state["index"]
                    copies[state["index"]] = 1
                    keep.append(state)
                else:
                    # Past the cut a family still needs one entry to absorb its own
                    # copies, or they would be counted against a state they are not
                    # copies of; it is simply never written out.
                    representative[key] = None
                    copies.setdefault(None, 0)
        # One trajectory per level, the kept states laid end to end: a viewer opens
        # `engaged_07.xyz` and has that level's whole shortlist as one sequence,
        # rather than twenty files to load and keep in order by hand.  The states run
        # in rank order, so the file reads best to worst and a frame's number divided
        # by the frames per state gives its rank.
        levelPath = os.path.join(outputDirectory, "engaged_%02d.xyz" % level)
        levelTrajectory = open(levelPath, "w")
        time = 0

        levelDirectory = outputDirectory if arguments.flat \
            else os.path.join(outputDirectory, "engaged_%02d" % level)
        if arguments.per_state:
            os.makedirs(levelDirectory, exist_ok=True)

        for rank, state in enumerate(keep, start=1):
            frames, problem = collectState(runDirectory, state)
            if problem is not None:
                incomplete.append((state["index"], problem))
                continue
            if len(frames) < 3:
                incomplete.append((state["index"], "no relaxed structure"))
            elif len(frames) < 4:
                incomplete.append((state["index"], "no freely relaxed structure"))

            firstFrame = time
            time = appendFrames(levelTrajectory, frames, state, rank,
                                arguments.sort_by, time)
            location = "%s:%d-%d" % (os.path.basename(levelPath),
                                     firstFrame, time - 1)

            if arguments.per_state:
                statePath = os.path.join(
                    levelDirectory, "rank%02d_state_%s.xyz" % (rank, state["index"]))
                with open(statePath, "w") as trajectory:
                    appendFrames(trajectory, frames, state, rank,
                                 arguments.sort_by, 0)
                location += "  " + os.path.relpath(statePath, outputDirectory)

            written += 1
            summary.write("%d  %d  %s  %s  %d  %s\n" % (
                level, rank, state["index"],
                "  ".join(str(state[column]) for column in passedThrough),
                copies.get(state["index"], 1),
                location))

        levelTrajectory.close()
        levelFiles.append((level, levelPath, time))

        best = keep[0][arguments.sort_by] if keep else float("nan")
        worst = keep[-1][arguments.sort_by] if keep else float("nan")
        collapsed = ""
        if dedup:
            merged = sum(count - 1 for index, count in copies.items()
                         if index is not None)
            collapsed = ",  %d copy(s) folded in" % merged
        print("  %2d node(s): %8d state(s) -> kept %2d,  %s %.6f .. %.6f%s"
              % (level, len(ranked), len(keep), arguments.sort_by, best, worst,
                 collapsed))

    summary.close()
    print("\nwrote %d state(s) into %d level trajectory(s) in %s"
          % (written, len(levelFiles), os.path.abspath(outputDirectory)))
    for level, levelPath, frames in levelFiles:
        print("  %s  %4d frame(s)" % (os.path.basename(levelPath), frames))
    print("summary: %s" % summaryPath)
    if incomplete:
        print("\n%d state(s) came out short:" % len(incomplete))
        for index, problem in incomplete[:20]:
            print("  state %s: %s" % (index, problem))
        if len(incomplete) > 20:
            print("  ... and %d more" % (len(incomplete) - 20))
    return 0


if __name__ == "__main__":
    sys.exit(main())
