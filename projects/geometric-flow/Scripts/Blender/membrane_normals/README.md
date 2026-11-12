# Membrane Normals

Blender add-on that puts arrows on a membrane mesh for rendering:

* the plain vertex normals (all the same length, one color),
* the mean curvature normals (same length, colored by the magnitude of the mean curvature),
* the Gaussian curvature and the smaller and larger principal curvature (arrows along the normals, colored by the
  signed value with a diverging colormap),
* the forces the simulation dumped for every energy term (colored by the force magnitude, optionally sized by it).

The arrows live in their own collection (`Normals - <mesh>` or `Forces - <mesh>`) and every parameter stays editable afterwards.

## Install

Blender 5.2 (also written for 4.2+, only tested on 5.2.2):

    cd projects/geometric-flow/Scripts/Blender
    blender --command extension build --source-dir membrane_normals --output-dir .

then Preferences > Get Extensions > (v) > Install from Disk > `membrane_normals-0.3.0.zip`.

0.2.0 changed how arrow sets are stored (direction and magnitude per point). Arrow sets made with 0.1.0 do not
work with it: delete them and add them again. Sets made with 0.2.x keep working in 0.3.0 (they keep their own node
group), without the newer fields.

## Use

1. Select the membrane mesh(es), 3D viewport sidebar > **Normals**, then **Normals**, **Mean Curvature**,
   **Gaussian**, **Min k** or **Max k** (the smaller and larger principal curvature).
   An arrow object appears in `Normals - <mesh>` and becomes active. Pressing the same button again only refreshes
   the existing set.
2. With the arrow object active the panel has:
   * **Offset**: distance of the arrow base from the mesh along the vertex normal. 0 (default) touches the mesh.
   * **Length**, **Thickness** (shaft radius / length), **Head Length** (fraction of the length). Every arrow has the same length.
     The default length is 1.5 mean edge lengths.
   * **Flip**: reverse the arrows (the offset still goes along the normal).
   * **Size by Magnitude** (colored sets): scale every arrow by its absolute value relative to the larger end of the
     color range, so **Length** is the arrow at that end and a zero force has no arrow.
   * **Density** (+ **Seed**): draw arrows on a random fraction of the vertices.
   * **Color by |H|, |F|, K, k min or k max** (colored sets): **Min / Max** is the range mapped onto the colormap
     (values outside are clamped). **Fit Range** sets 0 to the 98th percentile of the magnitudes (for K and the
     principal curvatures minus that value to plus it, so zero is the middle of the colormap), so a few spiky
     vertices do not wash out the rest. The **colormap** menu offers Viridis, Plasma, Inferno, Turbo, Cool to Warm.
     |H| and |F| start on Viridis, K and the principal curvatures on Cool to Warm; the two kinds have their own
     material (`Magnitude Arrows`, `Signed Arrows`), so a colormap you pick changes the sets that share that
     material, not the other kind. The panel also shows the value range found on the mesh.
   * **Hide Below** (colored sets): arrows whose absolute value is under this are not drawn, so a force of zero
     (a bead that does not touch the membrane) leaves no arrow. It starts at 1e-3 of the top of the color range
     for |H| and |F|, and at 0 for K and the principal curvatures (a curvature of zero is a result, not a missing
     force: a cylinder has K = 0 everywhere). Raise it to drop more, set it to 0 to draw everything. The panel says
     how many arrows it hides.
   * **Update from Mesh** measures the mesh again (see Limits); on a force set it is **Reload File**.
     The trash button deletes the arrow set.
3. Colors show in Material Preview, Rendered and the renders. Set the scene's View Transform to **Standard**
   for the exact colormap colors (AgX and Filmic shift them; the panel warns when it is not Standard).

## Forces from the simulation

1. In the config file of the run add `"dump_forces": true`. Every time `main_cluster` saves a mesh
   (`save_interval`, and the final state) it also writes, for every energy term, the force on every vertex:

       <run folder>/membrane_120.obj
       <run folder>/Forces/Bending_tan_120.txt
       <run folder>/Forces/Surface_tension_120.txt
       <run folder>/Forces/Volume_constraint_120.txt
       <run folder>/Forces/Bead_0_120.txt        (one per bead: the force of that bead on the membrane)

   The name is the term as written in the `Output_data.txt` header and the step is the one in `membrane_<step>.obj`.
   Each file has a `#` header line and then one `fx fy fz` line per vertex, in the vertex order of that `.obj`.
   The force is minus the energy gradient, the quantity the integrators add to the positions (`x += dt * F`).
   A term that is switched off (a constraint with constant 0) gives a file of zeros. The files are written by a
   separate evaluation of the forces at the saved state, so a run with `dump_forces` is bit for bit the same as
   without it, only slower at the saved steps and with larger folders (about 3 numbers per vertex per term per saved step).
2. In Blender show the mesh of that step: with an OBJ sequence set the frame to the one that shows `membrane_<step>.obj`
   (or import that file), and select the mesh.
3. Sidebar > Normals > **Load Forces**. The file browser opens; pick the exact file (or several with shift or ctrl),
   or pick one file and tick **All Terms of This Step** in the options to load every `<term>_<step>.txt` of that step
   in the same folder. Each file becomes an arrow set named `<mesh> <term> <step>` in `Forces - <mesh>`, so you
   can show, hide or render each term on its own and give each its own range.

The loader checks that the file has one line per vertex of the mesh, and, when the `membrane_<step>.obj` of the
run is next to the `Forces` folder, that the mesh in the scene is that state (a warning if it is not, for example
because the sequence is on another frame). The arrows use the base mesh of the object, so they are on every
vertex of it, also where a Membrane Slicer cut hides the surface.

## Definitions

* **Normals**: Blender's vertex normals. Arrows point along them.
* **Mean curvature normal**: the area gradient of the vertex, the same quantity as `vertexNormalMeanCurvature`
  in `core/src/geometry.cpp`: `1/2 sum_j (cot a_ij + cot b_ij) (x_i - x_j)`. It points the way the area grows
  (outward on a convex region, whatever the face orientation). The arrow is drawn along the vertex normal, flipped
  where this vector points against it. The magnitude is divided by twice the barycentric dual area, so it is the
  mean curvature H = (k1 + k2) / 2 and not the integrated vector: a sphere of radius R reads 1/R at every vertex
  (checked: R = 1 gives 1.0006 on average, a cylinder of radius 1 gives 0.5000, a plane 0).
  On a boundary (an open sheet or a slicer cut) the half ring area is used.
* **Gaussian curvature K**: angle defect over the barycentric dual area, as `angleDefect` in `core/src/geometry.cpp` and
  geometry-central have it (2 pi minus the angles around the vertex; pi minus them on a boundary, with the half ring
  area). 1/R^2 on a sphere. Independent of the face orientation, positive where the surface curves the same way in
  both directions (the rim of a disc, a sphere) and negative on saddles (the dimple of a biconcave disc, the inside
  of a torus).
* **Principal curvatures**: `H -/+ sqrt(H^2 - K)` with the H and K above, as the simulation's `principalCurvatures` and
  geometry-central's `vertexMin/MaxPrincipalCurvatures`; Min k is the smaller, Max k the larger. They have the sign
  of H, so with outward normals a sphere reads +1/R for both. Where the discrete H^2 is below K (a round patch)
  the root is 0 and both read H. Checked on a torus of R = 2, r = 0.5: K, k min and H are within 1e-3 of the exact
  values everywhere (including the sign change between the outer and inner equator), the angle defects sum to
  2 pi times the Euler characteristic (0 for the torus), and an open cylinder reads K = 0, k min = 0, k max = 1
  on its rims as well. The arrows are along the vertex normals (not the principal directions) and the color
  carries the value.
* **Force arrows**: along the force of the term on that vertex (the full 3D vector, not only its normal part),
  colored by its length. A vertex with no force has no direction: its arrow follows the normal.
* Units are those of the mesh in object space, so scaling the object does not change the readings.

## How it works

The arrow object's mesh is just one point per membrane vertex, with the attributes `ne_normal` (offset direction),
`ne_direction` (unit arrow direction) and `ne_magnitude` written by the add-on. The `Normal Arrows` Geometry Nodes
modifier turns the points into arrows (offset, thin out, shaft + head instanced, aligned and optionally scaled,
realized) and stores the normalized value as `ne_value` (`ne_magnitude` is signed for the curvatures), which the
`Magnitude Arrows` or `Signed Arrows` material feeds through a color ramp (edit it in the Shader editor, or use the
colormap menu; it is shared by the sets of the same material). Plain normals use the `Normal Arrows` material. The arrow object is parented to the membrane (identity parent inverse), so it
follows when you move the membrane.

Normal and mean curvature sets read the mesh as evaluated, so a Membrane Slicer cut (or any other modifier) is
respected: arrows only appear on what is left. Pass `use_modifiers=False` in the redo panel of **Add** to use the
base mesh instead.

## Limits

* The arrows are a snapshot: after changing the frame of an OBJ sequence (Stop-Motion-OBJ), editing the mesh or
  moving the slicer, press **Update from Mesh** (works on every selected arrow set). Force sets keep the file they
  were loaded from: to see another step, load that step's file.
* Vertex arrows only (no face normals), and one arrow shape (shaft + cone head, 12 sides, flat shaded).
* K and the principal curvatures are estimated per vertex from the triangles, so they scatter on the irregular
  meshes the simulation remeshes to: on the 453-vertex `Simplesphere.obj` (exact K = 1) the middle 96% of the
  vertices read K between 0.78 and 1.34 and the larger principal curvature goes up to 2.4 (exact value 1), the same
  as geometry-central gives. The root in `H +/- sqrt(H^2 - K)` amplifies the scatter where the surface is nearly
  round. The default color range (98th percentile of |value|) keeps the few worst vertices from washing out the
  rest; finer or more regular meshes are cleaner (the torus above is within 1e-3).
* Force files only include the energy terms of the simulation. The Lagrange multiplier forces of the constrained
  Newton integrators are not terms and are not written.
* The materials are plain Principled BSDFs, so the colors also respond to your lights. For flat, exact colors
  swap the shader for an Emission one in the Shader editor.
