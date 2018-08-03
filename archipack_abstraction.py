# -*- coding:utf-8 -*-

# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110- 1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

# ----------------------------------------------------------
# Author: Stephen Leger (s-leger)
#
# ----------------------------------------------------------
# noinspection PyUnresolvedReferences
import bpy
import time
import logging
logger = logging.getLogger("archipack")

"""
 Map between object classes name and layer index in 2.7 series
"""
layer_map = {
    "All":0,
    "reference_point":10,
    "wall":11,
    "custom_wall":11,
    "wall2":11,
    "slab":12,
    "slab_cutter":12,
    "floor":13,
    "floor_cutter":13,
    "molding":13,
    "fence":14,
    "stair":14,
    "truss":14,
    "kitchen":14,
    "kitchen_module":14,
    "handle":15,
    "window":15,
    "window_panel":15,
    "window_shutter":15,
    "blind":15,
    "door":15,
    "door_panel":15,
    "custom":15,
    "custom_part":15,
    "roof":16,
    "roof_cutter":16,
    "profile":17,
    "dimension_auto":18,
    "section":18,
    "section_target":18,
    "layout":18,
    "hole":19,
    "hybridhole":19,
    "custom_hole":19
}
"""
 Layer (layer Management) / collection names by class
"""
layer_names = {
    "All":"All",
    "reference_point":"Reference",
    "wall":"Walls",
    "custom_wall":"Walls",
    "wall2":"Walls",
    "slab":"Slabs",
    "floor":"Floors",
    "molding":"Floors",
    "fence":"Furnitures",
    "stair":"Furnitures",
    "truss":"Furnitures",
    "kitchen":"Furnitures",
    "kitchen_module":"Furnitures",
    "handle":"Openings",
    "window":"Openings",
    "window_panel":"Openings",
    "window_shutter":"Openings",
    "blind":"Openings",
    "door":"Openings",
    "door_panel":"Openings",
    "custom":"Openings",
    "custom_part":"Openings",
    "roof":"Roofs",
    "slab_cutter":"Slabs",
    "floor_cutter":"Floors",
    "roof_cutter":"Roofs",
    "dimension_auto":"2d",
    "section":"2d",
    "section_target":"2d",
    "layout":"2d",
    "profile":"Profile",
    "hole":"Holes",
    "hybridhole":"Holes",
    "custom_hole":"Holes"
}


class ArchipackLayerManager_27():
    
    """
     Layer manager for use in 2.7 series
    """
    def layer_by_name(self, context, layer_name):
        if layer_name in layer_map:
        
            return layer_map[layer_name]
        return -1
        
    def get_layer_index(self, context, o, layer_name):
        if layer_name is not None:
            index = self.layer_by_name(context, layer_name)
            
        else:    
            if o.data:
                for key in o.data.keys():
                    if "archipack_" in key:
                        layer_name = key[10:]
                        break
                        
            if layer_name is None:
                for key in o.keys():
                    if "archipack_" in key:
                        layer_name = key[10:]
                        break
            
            index = self.layer_by_name(context, layer_name)
                        
        if index > -1 and hasattr(context.scene, "namedlayers"):
            # setup layer name of Layer Management addon
            layer = context.scene.namedlayers.layers[index]
            if index < 10:
                name = "Layer0{}".format(index + 1)
            else: 
                name = "Layer{}".format(index + 1)
            if layer.name == name:
                layer.name = layer_names[layer_name]
            
        return index
        
    def add_to_layer(self, context, o, layer_name=None, default_layer=False):
        """
         Add object to layer using object class name
         layer_name: string use this layer_name instead of class one
         default_layer: add to default layer for all entity
         Dont add into default layer as select by layer using layer manager wont work
         This does trigger a DAG zero... report
        """
        index = self.get_layer_index(context, o, layer_name)
        tim = time.time()
        layers = {}
        if 20 > index > -1:
            layers[index] = True
            # ensure layer is visible
            if not context.scene.layers[index]:
                # does trigger DAG zero...
                context.scene.layers[index] = True
        
        if default_layer:
            index = self.layer_by_name(context, "All")
            if 20 > index > -1:
                layers[index] = True
                if not context.scene.layers[index]:
                    # does trigger DAG zero...
                    context.scene.layers[index] = True
                
        put_on_layers = lambda x: tuple((i in x) for i in range(20))
        o.layers[:] = put_on_layers(layers)        
        logger.debug("add_to_layer() :%.4f seconds", time.time() - tim)
        
    
class ArchipackLayerManager_28():
    
    """
     Layer manager for use in 2.8 series
    """
    def collection_by_name(self, context, collection_name):
        if collection_name in layer_map:
            return collection_name.capitalize()
        return None
        
    def add_to_layer(self, context, o, layer_name=None, default_layer=False):
        """
         Add object to layer using object class name
         layer_name: string use this layer_name instead of class one
         default_layer: add to default layer for all entity
        """
        if layer_name is not None:
            collection_name = self.collection_by_name(context, layer_name)
            
        else:    
            collection_name = None
            if o.data:
                for key in o.data.keys():
                    if "archipack_" in key:
                        collection_name = self.collection_by_name(context, key[10:])
                        break
            else:
                for key in o:
                    if "archipack_" in key:
                        collection_name = self.collection_by_name(context, key[10:])
                        break
        # @TODO:
        # add object to collection
        if 20 > index > -1:
            o.layers[index] = True
        
        if default_layer:
            collection_name = self.collection_by_name(context, "All")
            if 20 > index > -1:
                o.layers[index] = True
    
       
class ArchipackObjectsManagerBase():

    def _cleanup_datablock(self, d, typ):
        if d and d.users < 1:
            if typ == 'MESH':
                bpy.data.meshes.remove(d)
            elif typ == 'CURVE':
                bpy.data.curves.remove(d)
            elif typ == 'LAMP':
                bpy.data.lamps.remove(d)

    def _delete_object(self, context, o):
        d = o.data
        typ = o.type
        self.unlink_object_from_scene(context, o)
        bpy.data.objects.remove(o)
        self._cleanup_datablock(d, typ)

    def _delete_childs(self, context, o):
        for child in o.children:
            self._delete_childs(context, child)
        self._delete_object(context, o)

    def delete_object(self, context, o):
        """
          Recursively delete object and childs
          Cleanup datablock when needed
          @o: object to delete
        """
        if o is not None:
            self._delete_childs(context, o)

    def _duplicate_object(self, context, o, linked):
        new_o = o.copy()
        if o.data:
            if linked:
                new_o.data = o.data
            else:
                new_o.data = o.data.copy()
        self.link_object_to_scene(context, new_o)
        return new_o

    def _duplicate_childs(self, context, o, linked):
        p = self._duplicate_object(context, o, linked)
        for child in o.children:
            c = self._duplicate_childs(context, child, linked)
            c.parent = p
            # c.location = child.location.copy()
            c.matrix_local = child.matrix_local.copy()
            c.matrix_parent_inverse = child.matrix_parent_inverse.copy()
        return p

    def duplicate_object(self, context, o, linked):
        """
          Recursively duplicate object and childs
          @o: object to duplicate
          @linked : boolean linked duplicate
          return parent on success
        """
        if o is not None:
            return self._duplicate_childs(context, o, linked)
        return None

    def _link_object(self, src, o):
        if src.data:
            d = o.data
            typ = o.type
            o.data = src.data
            self._cleanup_datablock(d, typ)

    def _link_child(self, src, o):
        self._link_object(src, o)
        if len(src.children) == len(o.children):
            for i, child in enumerate(src.children):
                self._link_child(child, o.children[i])

    def link_object(self, src, o):
        """
         Recursievely link datablock
         @src: object source
         @o: object destination
         src and o parent child relationship must match
        """
        if src is not None:
            self._link_child(src, o)

    def get_topmost_parent(self, o):
        if o.parent:
            return self.get_topmost_parent(o.parent)
        else:
            return o
    
    def get_reference_point(self, o):
        if o.parent:
            return self.get_reference_point(o.parent)
        elif "archipack_reference_point" in o:
            return o
        else:
            return None
    
    def filter_by_class_name(self, o, cls):
        return o.data and cls in o.data
    
    def _get_objects_by_class_name(self, o, cls, res):
        if self.filter_by_class_name(o, cls):
            res.append(o)
        for c in o.children:
            self._get_objects_by_class_name(c, cls, res)
        
    def get_objects_by_class_name(self, o, cls):
        res = []
        ref = self.get_reference_point(o)
        if ref is not None:
            self._get_objects_by_class_name(ref, cls, res)
        return res
      
      
class ArchipackObjectsManager_27(ArchipackLayerManager_27, ArchipackObjectsManagerBase):
    """
      Provide objects and datablock utility
      Support meshes curves and lamps
      - recursive delete objects and datablocks
      - recursive clone linked
      - recursive copy

      Provide abstraction layer for blender 2.7 series
    """
    def scene_ray_cast(self, context, orig, vec):
        return context.scene.ray_cast(
            orig,
            vec)

    def get_scene_object(self, context, name):
        return context.scene.objects.get(name)

    def get_scene_objects(self, context):
        return context.scene.objects

    def select_object(self, context, o, active=False):
        """
         Select object and optionnaly make active
        """
        if o is not None:
            o.select = True
            if active:
                context.scene.objects.active = o
            
    def unselect_object(self, o):
        o.select = False
        
    def link_object_to_scene(self, context, o, layer_name=None, default_layer=False):
        tim = time.time()
        context.scene.objects.link(o)
        logger.debug("context.scene.objects.link(o) :%.4f seconds", time.time() - tim)
        self.add_to_layer(context, o, layer_name, default_layer)
        
    def unlink_object_from_scene(self, context, o):
        context.scene.objects.unlink(o)

    def is_visible(self, o):
        return o.hide is False

    def is_selected(self, o):
        return o.select

    def hide_object(self, o):
        o.hide = True

    def show_object(self, o):
        o.hide = False


class ArchipackObjectsManager_28(ArchipackLayerManager_28, ArchipackObjectsManagerBase):
    """
      Provide objects and datablock utility
      Support meshes curves and lamps
      - recursive delete objects and datablocks
      - recursive clone linked
      - recursive copy

     Provide abstraction layer for blender 2.8 series

    """
    def scene_ray_cast(self, context, orig, vec):
        return context.scene.ray_cast(
            context.scene.view_layers.active,
            "",
            orig,
            vec)

    def get_scene_object(self, context, name):
        return context.scene_collection.objects.get(name)

    def get_scene_objects(self, context):
        return context.scene_collection.objects

    def select_object(self, context, o, active=False):
        """
         Select object and optionnaly make active
        """
        if o is not None:
            act = context.object
            o.select_set(action='SELECT')
            if act is not None and not active:
                # reselect active object so it remains active one
                act.select_set(action='SELECT')

    def unselect_object(self, o):
        o.select_set(action='DESELECT')

    def link_object_to_scene(self, context, o, layer_name=None, default_layer=False):
        context.scene_collection.objects.link(o)
        self.add_to_layer(context, o, layer_name, default_layer)
        
    def unlink_object_from_scene(self, context, o):
        context.scene_collection.objects.unlink(o)

    def is_visible(self, o):
        return o.visible_get()

    def is_selected(self, o):
        return o.select_get()

    def hide_object(self, o):
        o.hide_viewport = True
        
    def show_object(self, o):
        o.hide_viewport = False
        