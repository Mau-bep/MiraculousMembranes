# Membrane Slicer

Blender add-on that cuts away a plane or a cake slice of a membrane, for every
frame of an OBJ sequence (made for the Stop-Motion-OBJ / OBJSequence add-on).
The cut is a modifier plus a normal object, so you can move, rotate and
keyframe it yourself.

## Install

Blender 5.2 (also written for 4.2+, only tested on 5.2.2):

    cd projects/geometric-flow/Scripts/Blender
    blender --command extension build --source-dir membrane_slicer --output-dir .

then Preferences > Get Extensions > (v) > Install from Disk > `membrane_slicer-0.1.1.zip`.

## Use

1. Select the membrane object(s) (the sequence object), 3D viewport sidebar > **Slicer** > **Create Slicer**.
   A wireframe wedge appears, opening toward the viewer. Its apex line is the slicer's local Z axis.
2. Move / rotate / scale the slicer object like anything else.
3. With the slicer selected, the panel has:
   * **Angle**: opening of the slice. **Plane (180°)** turns it into a single plane cut.
   * **Radius / Height**: how far the wedge reaches (**Fit to Targets** resizes it).
   * **Enabled**, **Limit to Frames** + Start / End: cut only a range of scene frames.
   * **Cut Membranes**: the objects being cut. Select more objects and press **Add Selected**.
   * **Keep Cut Faces**: keep the wall faces the boolean creates (off = open surface).
   * **Hole Tolerant** (open sheets) and **Self Intersections** (membranes that overlap themselves, slower).
4. Everything (slicer transform, Angle, Enabled, frame range) can be keyframed with `I` over the field.

## How it works

No frame handlers and no changes to the sequence's mesh data. Two Geometry Nodes groups:

* `MS Cutter` (on the slicer): builds the closed wedge solid and stores the settings on it.
* `MS Cut` (first modifier on each membrane): reads the slicer's geometry and settings, subtracts the
  wedge with an exact mesh boolean, deletes the cap faces and triangulates. Because it is a modifier,
  it follows whatever mesh the sequence addon puts on the object each frame.

The slicer is hidden from renders. Delete it with **Delete Slicer** (also removes the modifiers).

## Limits

* Cost is an exact boolean per frame: roughly 30 to 300 ms for 1 to 4k triangle membranes.
* Only objects you add are cut (a bead in the same scene stays whole).
* Output is triangles; faces split by the cut can be thin slivers at the rim.
