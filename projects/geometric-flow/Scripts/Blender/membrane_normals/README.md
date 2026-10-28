# Membrane Normals

Blender add-on that puts arrows on a membrane mesh for rendering:

* the plain vertex normals (all the same length, one color),
* the mean curvature normals (same length, colored by the magnitude of the mean curvature),
* the forces the simulation dumped for every energy term (colored by the force magnitude, optionally sized by it).

The arrows live in their own collection (`Normals - <mesh>` or `Forces - <mesh>`) and every parameter stays editable afterwards.

## Install

Blender 5.2 (also written for 4.2+, only tested on 5.2.2):

    cd projects/geometric-flow/Scripts/Blender
    blender --command extension build --source-dir membrane_normals --output-dir .

then Preferences > Get Extensions > (v) > Install from Disk > `membrane_normals-0.2.0.zip`.

0.2.0 changed how arrow sets are stored (direction and magnitude per point). Arrow sets made with 0.1.0 do not
work with it: delete them and add them again.

## Use

1. Select the membrane mesh(es), 3D viewport sidebar > **Normals**, then **Normals** or **Mean Curvature**.
   An arrow object appears in `Normals - <mesh>` and becomes active. Pressing the same button again only refreshes
   the existing set.
2. With the arrow object active the panel has:
   * **Offset**: distance of the arrow base from the mesh along the vertex normal. 0 (default) touches the mesh.
   * **Length**, **Thickness** (shaft radius / length), **Head Length** (fraction of the length). Every arrow has the same length.
     The default length is 1.5 mean edge lengths.
   * **Flip**: reverse the arrows (the offset still goes along the normal).
   * **Size by Magnitude** (colored sets): scale every arrow by its magnitude on the color range, so **Length** is the
     arrow at Range Max and a zero force has no arrow.
   * **Density** (+ **Seed**): draw arrows on a random fraction of the vertices.
   * **Color by |H|** or **|F|** (colored sets): **Min / Max** is the magnitude range mapped onto the colormap
     (values outside are clamped). **Fit Range** sets 0 to the 98th percentile of the magnitudes, so a few spiky
     vertices do not wash out the rest. The **colormap** menu offers Viridis, Plasma, Inferno, Turbo, Cool to Warm.
     The panel also shows the magnitude range found on the mesh.
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
* **Force arrows**: along the force of the term on that vertex (the full 3D vector, not only its normal part),
  colored by its length. A vertex with no force has no direction: its arrow follows the normal.
* Units are those of the mesh in object space, so scaling the object does not change the readings.

## How it works

The arrow object's mesh is just one point per membrane vertex, with the attributes `ne_normal` (offset direction),
`ne_direction` (unit arrow direction) and `ne_magnitude` written by the add-on. The `Normal Arrows` Geometry Nodes
modifier turns the points into arrows (offset, thin out, shaft + head instanced, aligned and optionally scaled,
realized) and stores the normalized magnitude as `ne_value`, which the `Magnitude Arrows` material feeds through a
color ramp (edit it in the Shader editor, or use the colormap menu; it is shared by all colored sets). Plain normals
use the `Normal Arrows` material. The arrow object is parented to the membrane (identity parent inverse), so it
follows when you move the membrane.

Normal and mean curvature sets read the mesh as evaluated, so a Membrane Slicer cut (or any other modifier) is
respected: arrows only appear on what is left. Pass `use_modifiers=False` in the redo panel of **Add** to use the
base mesh instead.

## Limits

* The arrows are a snapshot: after changing the frame of an OBJ sequence (Stop-Motion-OBJ), editing the mesh or
  moving the slicer, press **Update from Mesh** (works on every selected arrow set). Force sets keep the file they
  were loaded from: to see another step, load that step's file.
* Vertex arrows only (no face normals), and one arrow shape (shaft + cone head, 12 sides, flat shaded).
* Force files only include the energy terms of the simulation. The Lagrange multiplier forces of the constrained
  Newton integrators are not terms and are not written.
* The materials are plain Principled BSDFs, so the colors also respond to your lights. For flat, exact colors
  swap the shader for an Emission one in the Shader editor.
