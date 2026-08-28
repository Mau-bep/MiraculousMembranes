bl_info = {
    "name": "My Addon",
    "author": "Mau",
    "version": (1, 0),
    "blender": (2, 80, 0),
    "location": "View3D > Tool",
    "description": "We will render the mesh eventually",
    "warning": "",
    "doc_url": "",
    "category": "Add Mesh",
}



import bpy
import math

def cylinder_between(x1, y1, z1, x2, y2, z2, r):

  dx = x2 - x1
  dy = y2 - y1
  dz = z2 - z1    
  dist = math.sqrt(dx**2 + dy**2 + dz**2)

  bpy.ops.mesh.primitive_cylinder_add(
      radius = r, 
      depth = dist,
      location = (dx/2 + x1, dy/2 + y1, dz/2 + z1)   
  ) 

  phi = math.atan2(dy, dx) 
  theta = math.acos(dz/dist) 

  bpy.context.object.rotation_euler[1] = theta 
  bpy.context.object.rotation_euler[2] = phi 
  
#  
curves = bpy.data.curves
objects = bpy.data.objects
scene = bpy.context.scene

obj_name = "EdgesObject"
obj = objects[obj_name]

def make_tubes(obj, bevel_depth=0.0008, resolution=3):

    mesh = obj.data
    curve_name = 'TubesCurve'
    mat = obj.matrix_world
    # if exists, pick up else generate a new one
    cu = curves.get(curve_name, curves.new(name=curve_name, type='CURVE'))
    cu.dimensions = '3D'
    cu.fill_mode = 'FULL'
    cu.bevel_depth = bevel_depth
    cu.bevel_resolution = resolution
    cu_obj = objects.get(curve_name, objects.new(curve_name, cu))

    # break down existing splines entirely.
    if cu.splines:
        cu.splines.clear()

    # and rebuild
    verts = mesh.vertices
    for e in mesh.edges:
        idx_v1, idx_v2 = e.vertices
        v0, v1 = mat @ verts[idx_v1].co, mat @ verts[idx_v2].co
        full_flat = [v0[0], v0[1], v0[2], 0.0, v1[0], v1[1], v1[2], 0.0]

        # each spline has a default first coordinate but we need two.
        segment = cu.splines.new('POLY')
        segment.points.add(1)
        segment.points.foreach_set('co', full_flat)

    if not curve_name in scene.objects:
#        scene.objects.link(cu_obj)
        bpy.context.collection.objects.link(cu_obj)  


make_tubes(obj)


def point_cloud(ob_name, coords):
    """Create point cloud object based on given coordinates and name."""
    me = bpy.data.meshes.new(ob_name + "Mesh")
    ob = bpy.data.objects.new(ob_name, me)
    me.from_pydata(coords, [], [])
    ob.show_name = True
    me.update()
    return ob

def vertexcoords(obj):
    """Returns median center coordinates for each face of given mesh object."""
    if obj.type == 'MESH':
        import bmesh
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        return [obj.matrix_world @ v.co for v in obj.data.vertices]
    else:
        return [(0.0, 0.0, 0.0)]


ob = bpy.context.active_object
pc = point_cloud(ob.name + "-pointcloud", vertexcoords(ob))

# Link object to the active collection
bpy.context.collection.objects.link(pc)



class WM_OT_HelloWorld(bpy.types.Operator):
    bl_idname = "wm.hello_world"
    bl_label = "Minimal Operator"
    bl_options = {'REGISTER'}

    # Operator user properties, should be assigned using a single colon :
    # instead of using an equal sign = in Blender 2.8
    report_flag: bpy.props.BoolProperty(
        name = "Report",
        default = True)

    @classmethod # Will never run when poll returns false
    def poll(cls, context):
        return context.object

    def invoke(self, context, event): # Used for user interaction
        wm = context.window_manager
        return wm.invoke_props_dialog(self)

    def draw(self, context): # Draw options (typically displayed in the tool-bar)
        row = self.layout
        row.prop(self, "report_flag", text="Report Hello World")

    def execute(self, context): # Runs by default 
        if self.report_flag:
            self.report({'INFO'}, "Hello World")
        else:
            print ("Hello World")
            
            
        #We want to do somethign more
        obj = context.object
        data = obj.data
        vertices = data.vertices
        mat = obj.matrix_world
#        make_tubes(obj)
#        for e in data.edges:
#            idx_v1, idx_v2 = e.vertices
#            pos_v1 = mat @ vertices[idx_v1].co
#            pos_v2 = mat @ vertices[idx_v2].co 
#            cylinder_between(pos_v1[0],pos_v1[1],pos_v1[2],pos_v2[0],pos_v2[1],pos_v2[2],0.02)
        for v in data.vertices:
            pos = mat @ v.co
            bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=3, radius=0.003, enter_editmode=False, align='WORLD', location=(pos[0], pos[1], pos[2]), scale=(1, 1, 1))
        
        print(type(vertices[0].co))
#        addConeBetweenPoints(vertices[0].co,vertices[1].co,0.2,0.0)
        print(vertices[0].co)
        return {'FINISHED'}


class TestPanel(bpy.types.Panel):
        bl_label= "TestPanel "
        bl_idname = "PT_TestPanel"
        bl_space_type = "VIEW_3D"
        bl_region_type = 'UI'
        bl_category = 'Addon'
        
        
        def draw(self,context):
            
            layout = self.layout
            
            row = layout.row()
            row.label(text = "Sometext",icon="GHOST_ENABLED")
            row = layout.row()
            row.operator("mesh.primitive_cube_add")
            row = layout.row()
#            layout.operator.hello_world()

class PanelA(bpy.types.Panel):
        bl_label= "PanelA"
        bl_idname = "PT_PanelA "
        bl_space_type = "VIEW_3D"
        bl_region_type = 'UI'
        bl_category = 'Addon'
        bl_parent_id = "PT_TestPanel"
        bl_options ={"DEFAULT_CLOSED"}
        def draw(self,context):
            
            layout = self.layout
            obj = context.object
            row = layout.row()
            row.label(text = "Panel A FTW",icon="GHOST_ENABLED")
            row = layout.row()
            row.prop(obj,'scale')
#            vertices = obj.data.vertices
            row = layout.row()
            
#            bpy.ops.wm.hello_world('INVOKE_DEFAULT')
            row.operator("wm.hello_world")
#            plain_verts = [obj.matrix_world @ vert.co for vert in vertices]
#            print(plain_verts[0])
            
            
            
        

def register():
    bpy.utils.register_class(TestPanel)
    bpy.utils.register_class(PanelA)
    bpy.utils.register_class(WM_OT_HelloWorld)
    
def unregister():
    bpy.utils_unregister_class(TestPanel)
    bpy.utils.unregister_class(PanelA)
    bpy.utils.unregister_class(WM_OT_HelloWorld)
    
    
if __name__ == "__main__":
    register()
    
    
#bpy.ops.wm.hello_world('INVOKE_DEFAULT')
    
    
