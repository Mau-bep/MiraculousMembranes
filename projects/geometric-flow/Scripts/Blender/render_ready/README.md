# RenderReady

Blender add-on that dresses a mesh for a render: the mesh itself, a wireframe on top of it, and a point on every
vertex, each with a material you choose.

## Install

Blender 5.2 (also written for 4.2+, only tested on 5.2.2):

    cd projects/geometric-flow/Scripts/Blender
    blender --command extension build --source-dir render_ready --output-dir .

then Preferences > Get Extensions > (v) > Install from Disk > `render_ready-0.1.1.zip`.

## Use

1. Select the mesh (an OBJ-sequence object works), 3D viewport sidebar > **RenderReady** > **Make Render Ready**.
   Two copies appear in a collection `RenderReady - <mesh>`:
   * `<mesh> wire`: your mesh with a **Wireframe** modifier (Replace Original, so only the wires show).
   * `<mesh> points`: your mesh with a **Geometry Nodes** modifier (`RenderReady Points`) that does
     **Mesh to Points** on its vertices and then **Set Material**.
2. With the mesh (or either copy) selected, the panel has:
   * **Mesh Material**: material of the original mesh.
   * **Wire and Points Material**: one material for both the wireframe and the points.
     Both pickers have a **New** button to make a material on the spot.
   * **Wire Width**: thickness of the wireframe.
   * **Point Radius**: radius of the points.
   * **Even Thickness**: the Wireframe modifier option of the same name, off by default. On, it keeps the wire width
     constant at sharp corners, but on the thin sliver triangles that a Membrane Slicer cut leaves it throws wire
     vertices far off the surface, which show up as long stray lines across the picture (measured: up to 12 units on a
     membrane of radius 3). Leave it off for cut membranes.
   * **Remove Copies** deletes the two copies (the mesh keeps its Mesh Material).
3. Width and radius start at 0.06 and 0.1 of the mean edge length. Edit them, the materials and the mesh as you like:
   everything is live.

Points render as spheres in Cycles (tested). EEVEE has drawn point clouds since 4.2, but I have not tried it here.

## On top of an OBJ sequence (videos)

The copies **share the mesh data** of the original and nothing is applied: the Wireframe and the Geometry Nodes
modifiers are live. The OBJ-sequence add-on rewrites the one mesh in place on every frame change, so the wireframe
and the points follow every frame of the sequence, for stills and for videos. The same goes for anything else that
changes the mesh: edit mode, shape keys, and the modifiers the original has when you press Make Render Ready,
which the copies carry before their own wireframe or points modifier. So with a Membrane Slicer already on the mesh
the wireframe and the points are cut the same way (tested). If you add the slicer after RenderReady, select the
two copies together with the slicer and press **Add Selected** in the slicer panel.

The copies are parented to the original with no offset, so they move with it. Because the OBJ-sequence add-on would
also update every object that carries its sequence settings, RenderReady switches those settings off on the two
copies (the original keeps them), so each frame is read once.

## Materials

The mesh data is shared, so the materials are set on the objects (object-linked slots), not on the mesh:
the original keeps its own slots untouched except for the Mesh Material override, and the wireframe gets its
material without changing the points or the original. A mesh with no material slot gets one empty slot (on the
shared mesh) so there is something to assign to. Clearing a material field goes back to what the mesh itself has.
The Mesh Material goes on every slot of the original.

## Limits

* Edit the settings in the RenderReady panel: the panel pushes them onto the modifiers, so a change made directly
  in a modifier is overwritten the next time a setting changes.
* The points copy shows only points, the wire copy only wires: to see the surface keep the original visible.
* A mesh without edges takes 0.1 as its edge length for the first width and radius.
