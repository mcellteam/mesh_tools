# meshmorph_patches

Local patches against **mesh_tools**, needed to build and run `meshmorph` /
`meshalyzer` on this project's ROI meshes.

```
upstream   https://github.com/mcellteam/mesh_tools.git
base       205ebbf99a6211b7508c6a24f5cc5dcbbd825af1  (2025-03-11)
applied    2026-08-14, macOS / Apple clang / libc++, C++11
```

Three patches. Two are build fixes that any modern toolchain needs; one is a
performance fix that only matters at our scale, but at our scale it is the
difference between a run finishing and a run never finishing.

| patch | kind | what it touches |
| --- | --- | --- |
| 0001 | build | `meshmorph/container.cc`, `meshmorph/gain_schedule.cc`, `meshalyzer/meshalyzer.h` |
| 0002 | perf | `meshmorph/state.cc` |
| 0003 | perf | `meshmorph/intersecting_faces.h`, `meshmorph/state.cc`, `meshmorph/meshmorph.cc` |

## Applying

```sh
git clone https://github.com/mcellteam/mesh_tools.git
cd mesh_tools
git checkout 205ebbf
git am /path/to/meshmorph_patches/*.patch      # or: patch -p1 < each
cd meshmorph  && make -j8
cd ../meshalyzer && make -j8
```

### Applying to the `mcell4_dev` branch instead

These were authored against `205ebbf`, but all six files they touch are
**byte-identical** on `mcell4_dev` (7 commits ahead as of 2026-08-14), so
`git am` applies them cleanly there too. That is where they are installed
locally: branch `meshmorph-patched`, off
`/Users/bartol/src/mcell-git/mcell4_head/mesh_tools`.

### The OBJ converters DO build — use them, not a hand-rolled script

An earlier version of this file said `obj2mesh` and `mesh2obj` could not be
built. That was true of `205ebbf` and is **no longer true**: `mcell4_dev`
carries `a33533b`, `dc1d326`, `ce9dbbe` and `1cae62e`, which modernise
`obj2mesh`, `mesh2obj` and `meshscale` for current C, lex and bison. All three
build with no patching.

Our meshes are exported in microns and every meshmorph distance parameter is in
nanometres, so the conversion needs a scale:

```sh
obj2mesh in.obj 2>/dev/null | meshscale 1000 1000 1000 - > out.mesh
```

`obj2mesh` writes its `polygon mesh: N vertices & M polygons` banner to
**stderr**, so the `.mesh` on stdout is clean. `meshscale` reading `-` as stdin
comes from `d935565`. Verified 2026-08-14 against an independently produced
`.mesh`: identical face lists, maximum vertex difference 0.0 nm.

## 0001 — build: two C++11 errors

`container.cc:984` and `gain_schedule.cc:108` stream an `ofstream` into
`cout`; both meant to print the filename. Pre-C++11 this compiled through
`basic_ios`'s implicit `void*` conversion, which C++11's `explicit operator
bool` removed.

`meshalyzer.h` passes an explicit `std::allocator<V>` as the allocator template
argument to `std::map` / `multimap` / `unordered_map`. An associative
container's allocator `value_type` must be `pair<const K,V>`, and libc++
static-asserts on the mismatch. The defaults were always the correct types, so
the 21 explicit arguments are simply dropped.

## 0002 — perf: the per-move symmetry assert

`State::assignNewVertexCoords` runs once per vertex move and asserted
`Intersecting_Faces::intFacesAreSymmetric()`, which traverses the whole
intersecting-faces map doing a hash lookup per entry — **O(intersections) per
move**. The identical call inside `meshmorph.cc`'s move loop was already
commented out upstream, so this was the surviving copy.

The makefile does not define `NDEBUG`, so the assert is live in the `-O3`
build. It is invisible on models with few intersections and dominant on models
with many.

## 0003 — perf: constant-time emptiness tests

The one that mattered. `getCountOfIntFaces()` traverses the entire
intersecting-faces map and returns `count/2`. Two call sites used it purely to
compare against zero, both on the per-move path:

```
state.cc      if (getCountOfIntFaces(false) > 0)    every move
meshmorph.cc  if (getCountOfIntFaces(false) == 0)   every successful move,
                                                    live only under --no_int_prevention
```

So `--no_int_prevention` — the mode needed to separate interpenetrating objects
— paid **two full O(intersections) traversals per vertex move**. Replaced with
`hasNoIntersections()`, i.e. `intf.empty()`.

**Behavioural note, the only one in this patch set.** For the `meshmorph.cc`
self-tightening check, `empty()` is *stricter* than the old test: it requires
all intersections to be gone, including self-intersections, where
`getCountOfIntFaces(false)` counted only cross-object ones. It can therefore
only delay re-enabling strict face-intersection prevention, never enable it
early. That direction was chosen deliberately.

## Measurements

Same input, same flags, same machine. 809 objects, 2,589,343 vertices,
5,175,950 faces, 2,961,004 intersections; `--no_int_prevention -n 1 -g 10000
-K 1800`.

| build | moves | duration | rate |
| --- | --- | --- | --- |
| as shipped | 1,443 | 1,741 s | 0.83 /s |
| + 0002 | 1,872 | 1,764 s | 1.06 /s |
| + 0003 | 10,083 | 35 s | **286 /s** |
| rebuilt from these patches | 10,042 | 27.9 s | 360 /s |

0002 alone bought 28% because it removed one of three copies of the same O(N)
walk; 0003 removed the other two.

**Correctness check.** Startup state is bit-identical across every run —
1,502,800 nonnice vertices, 2,961,004 face intersections, minimum edge angle
0.4722981867 degrees — and intersections removed per move is unchanged
(~0.3/move). The patches change cost, not behaviour.

## Verifying after a re-apply

Run any `--measure_ecw_and_exit` or short `--no_int_prevention` job and confirm
the three startup numbers above. If they move, something other than cost has
changed.

Note `--measure_ecw_and_exit` **exits with status 1 on success** —
`meshmorph.cc` does `writeSepDistances(1); exit(1);`. Do not read that as a
failure.
