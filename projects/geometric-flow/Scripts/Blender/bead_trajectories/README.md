# Bead Trajectories

Blender add-on that turns the bead files of a simulation into animated objects: one sphere per bead that follows its
trajectory, and, if you want them, an arrow along the force on it. No OBJ sequence of spheres is needed.

## Install

Blender 5.2 (also written for 4.2+, only tested on 5.2.2):

    cd projects/geometric-flow/Scripts/Blender
    blender --command extension build --source-dir bead_trajectories --output-dir .

then Preferences > Get Extensions > (v) > Install from Disk > `bead_trajectories-0.1.2.zip`.

## Use

1. 3D viewport sidebar > **Beads** tab. Set **Results Folder** to the folder of one run (the one that holds
   `Bead_0_data.txt`, `Bead_1_data.txt`, ..., `Input_file.json` and the `membrane_<step>.obj` files).
2. **Start Frame**: the frame of the first mesh of the membrane OBJ sequence (1 by default).
3. **Load Beads**. A collection `Beads - <run>/<folder>` appears with, for every bead, an empty `Bead i` (keyframed
   location) and a child sphere `Bead i sphere` with its own material (the colours are spread over the hue circle; change
   them in the materials). All spheres share one mesh, scaled by the radius of the bead. **Sphere Detail** is the
   number of icosphere subdivisions (3 = 642 vertices). **Radius** 0 reads the radius of every bead from
   `Input_file.json`; a value overrides it.
4. **Same material for all beads** (shown once the beads exist): pick or create a material in the **Material** field and press
   **Apply to All Beads**; every sphere gets it (the arrows keep theirs). **Update Beads** brings back one colour per bead.
5. Pressing the Load button again (now **Update Beads**) rebuilds the collection from the files, so you can load a run that
   is still going and update it. **Remove Beads** deletes the collection.

## Frames

The frames are those of the saved meshes: the n-th `membrane_<step>.obj` (sorted by step, the order of an OBJ
sequence) is frame `Start Frame + n`, and the bead is placed at its row of that step in `Bead_<i>_data.txt`. If the bead
file has no row for a step the last earlier row is used, as when a run is resumed; if a step appears twice (a restarted
run) the last row wins. A folder without meshes takes the steps of the first bead file. Between two frames the
beads move in straight lines (linear keyframes).

## Force arrows

Tick **Force Arrows** (before Load Beads): every bead gets an arrow along the force on it (the `Fx Fy Fz` columns of
its file), starting at the surface of the sphere, in the material `Bead Force Arrows` (red). Its length is
**Arrow Length** bead radii for a typical force (the 95th percentile of all the forces of the run), proportional to
the force, and never longer than **Arrow Cap** times that. Forces below a thousandth of the typical one (a bead
that nothing touches) show no arrow. The arrow is a child of a pivot empty that carries the rotation keys and has its
own scale keys, so it renders in Cycles and EEVEE like any object. Changing Arrow Length or Cap needs **Update Beads**.
The arrows were checked for direction, base position and length on real data, and rendered once in Cycles;
they have not been looked at in a long animation.

## Moving the whole scene

The beads are keyframed, so moving a bead by hand is undone by its own animation. To turn or move the whole scene (beads,
membrane, slicer, ...) put everything under one pivot: select the membrane (the sequence object, plus the slicer or the camera
if you want them to go along; the RenderReady copies follow their membrane by themselves) and press **Group With Pivot**
in the **Move everything together** box. An empty `Scene Pivot` appears at the centre of everything and becomes the active
object; every selected object and every bead is parented to it without moving. Now **G**, **R** and **S** on the pivot move,
rotate and scale the whole scene, and it stays there for every frame, because the bead keyframes are relative to the pivot.
(Keep Auto-Keying off while you do it, or the pivot gets keyframes.) Pressing **Update Beads** later puts the new beads under
the same pivot. Objects that already have a parent of their own (other than the pivot) are left alone, so a RenderReady copy
is not taken out of its membrane. To undo it, delete the pivot with **Object > Delete Global**
after clearing the parents with **Alt+P > Clear and Keep Transformation**.

## On top of the membrane

The beads are ordinary objects, independent of the membrane sequence. With Stop-Motion-OBJ, set its start frame
to the same **Start Frame** and both play together; the Membrane Slicer cuts only the membrane, so a bead inside a cut
region is shown (move it by hand or hide the collection for the cut). RenderReady works on the membrane as before.

## Limits

* Radii come from `Input_file.json` of the run; beads with different radii work, the colours do not depend on them.
* The force is the total force on the bead as the simulation saved it at that step; steps without a saved row reuse the
  earlier one.
* Interpolation is linear in position, and the arrow rotates by interpolating its quaternion between frames.
