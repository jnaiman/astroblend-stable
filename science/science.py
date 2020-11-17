"""


Author: Jill Naiman <jnaiman@gmail.com>
Affiliation: UCSC
Homepage: www.astroblend.com
License:
  Copyright (C) 2012-2013 Jill Naiman.  All Rights Reserved.

  This file is part of blender science.

  science_blender is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

#!/usr/bin/python
# Filename: science.py
import bpy
import sys
from math import *
import os

from lighting import Lighting
from camera import Camera
from render import Render
from movie import Movies
import simpleobjects

from scienceutils import deselect_all, makeMaterial, setMaterial, delete_object, hide_object, unhide_object
import scienceutils

#import sys
#import cStringIO

# load files, pars what sort of file they are
class Load(object):
    def __init__(self, filelist, scale=(1.0,1.0,1.0), halo_sizes=[0.0008], particle_num=1, particle_colors=[(1,1,1)], 
                 isosurface_value = 1e-27, surf_type='sphere', radius = 10.0, 
                 radius_units = "mpc", surface_field="density",  
                 meshname = 'Allen', factor=1.0, 
                 transparency = 1.0, dist_fac = None,
                 color_field = None, emit_field = None, color_map = "algae", 
                 color_log = True, emit_log = True, plot_index = None, 
                 color_field_max = None, color_field_min = None, 
                 emit_field_max = None, emit_field_min = None, emissivity_units = None, 
                 n_ref=8, coords_tag = None):
        ifile = 0
        location = (0,0,0)
        particle_hide = []
        if isinstance(halo_sizes, tuple):
            hs = []
            for h in halo_sizes:
                hs.append(h)
            halo_sizes = hs
        # format halo sizes
        if isinstance(halo_sizes, float):
            halo_sizes = [halo_sizes]
        # reformat if a single string
        if isinstance(filelist,str):
            filelist = [filelist]
        # how long is the file?
        filelistlength = len(filelist)
        # now, figure out the file type based on the extention... for now
        if filelist[ifile].find('.obj') is not -1:
            filetype = 'obj'
        elif filelist[ifile].find('.txt') is not -1:
            filetype = 'sphtxt'
        else: # default behavior is to assume we are dealing with yt files
            filetype = 'yt' 
        # now, load if single file in all cases
        #  -> only file for single case, first file for many files
        if filetype is 'obj':
            name = self.load_obj(filelist[ifile])
            ds = []
        elif filetype is 'sphtxt':
            scale = [scale]
            ds = []
            if len(halo_sizes) < particle_num:
                hs = halo_sizes[0]
                halo_sizes = []
                for i in range(0,particle_num):
                    halo_sizes.append(hs)
            for i in range(0,particle_num):
                scale.append(scale[0])
            name = self.new_particle_set(filelist[ifile], particle_num, 
                                         particle_colors, halo_sizes, scale)      
        elif filetype is 'yt':
            # figure out what sorts of yt data we are dealing with - amr or sph?
            from yt import load
            ds = load(filelist[ifile])
            h = ds.h # we have to do this to get the dataset type to load... for real
            if (ds.dataset_type == 'flash') or (ds.dataset_type == 'enzo') or (ds.dataset_type == 'athena') or (ds.dataset_type == 'enzo_packed_3d'):
                filetype = 'ytamr' # this supports surfaces
                name = [meshname + '_0']
                ds = self.import_ytsurface(filelist[ifile],  isosurface_value, surf_type, radius, 
                         radius_units, surface_field, meshname, 
                         transparency, dist_fac, color_field, emit_field, color_map, 
                         color_log, emit_log, plot_index, color_field_max, color_field_min, 
                         emit_field_max, emit_field_min, emissivity_units)
                # futz with naming if multiple surfaces
                if isinstance(isosurface_value, float):
                    isosurface_value = [isosurface_value]
                    location = [location]
                    filelist = [filelist]
                    filetype = [filetype]
                    filelistlength = [filelistlength]
                    scale = [scale]
                if len(isosurface_value) > 1:
                    location = [location]
                    filelist = [filelist]
                    filetype = [filetype]
                    filelistlength = [filelistlength]
                    scale = [scale]
                    for i in range(1,len(isosurface_value)):
                        name.append(meshname + '_' + str(i))
                        location.append(location[0])
                        filelist.append(filelist[0])
                        filetype.append(filetype[0])
                        filelistlength.append(filelistlength[0])
                        scale.append(scale[0])
            elif (ds.dataset_type is 'tipsy') or (ds.dataset_type is 'gadget_hdf5'):
                filetype = 'ytsph'
                # tipsy default
                if ds.dataset_type is 'tipsy':
                    if coords_tag is None:
                        coords_tag = ('Gas', 'Coordinates')
                    if color_field is None:
                        color_field = ('Gas', 'Temperature')
                if ds.dataset_type is 'gadget_hdf5':
                    if coords_tag is None:
                        coords_tag = ('PartType0', 'Coordinates')
                    if color_field is None:
                        color_field = ('PartType0', 'InternalEnergy')
                name, ds, halo_sizes = self.import_ytsph(filelist[ifile], halo_sizes, color_field, 
                                                         color_map, color_log, n_ref, scale, coords_tag)
            else:
                print("SORRY!  Other yt data formats aren't supported yet! :(")
        else:
            print("DATA FORMAT NOT SUPPORTED!!!!")
        if isinstance(name, str):
            name = [name]
            location = [location]
            filelist = [filelist]
            filetype = [filetype]
            filelistlength = [filelistlength]
            scale = [scale]  
        # or if sph, for # of particles
        # please note: THIS IS SLIGHTLY HACKTACULAR
        if filetype is 'sphtxt':
            #if isinstance(filelist,str):
            filelist=[filelist]
            filelistlength = [filelistlength]
            filetype = [filetype]
            location = []
            particle_hide = []
            #scale = [scale]
            for i in range(0,particle_num):
                location.append((0,0,0))
                particle_hide.append(False)
                #scale.append(scale[0])
        if filetype is 'ytsph':
            filelist=[filelist]
            filelistlength = [filelistlength]
            filetype = [filetype]
            location = []
            particle_hide = []
            for i in range(0,len(name)):
                location.append((0,0,0))
                particle_hide.append(False)
                scale.append(scale[0])
        self.filetype = filetype
        self.__name = name # return the name of the thing too
        self.name = name
        self.location = location
        self.filelist = filelist
        self.filelistlength = filelistlength
        self.scale = scale
        self.__ifile = 0
        self.isosurface_value = isosurface_value
        self.surf_type = surf_type
        self.radius = radius
        self.radius_units = radius_units
        self.surface_field = surface_field
        self.transparency = transparency
        self.dist_fac = dist_fac
        self.color_field = color_field
        self.emit_field = emit_field
        self.color_map = color_map
        self.color_log = color_log
        self.emit_log = emit_log
        self.plot_index = plot_index
        self.color_field_max = color_field_max
        self.color_field_min = color_field_min
        self.emit_field_max = emit_field_max
        self.emit_field_min = emit_field_min
        self.emissivity_units = emissivity_units
        self.particle_num = particle_num
        self.particle_colors = particle_colors
        self.__halo_sizes = halo_sizes
        self.halo_sizes = halo_sizes 
        self.particle_hide = particle_hide
        self.ds = ds
        self.color_field = color_field
        self.color_map = color_map
        self.n_ref = n_ref
        self.coords_tag = coords_tag
        

    # for loading an objfile
    def load_obj(self, filename):
        bpy.ops.import_scene.obj(filepath=filename) # import an obj file
        #  first, find name that will be given to mesh
        xnn2 = 0
        xnn = 0
        # you only want the file name, not all the directory stuff
        while xnn != -1:
            xnn2 = xnn
            xnn = filename.find('/',xnn+1)
        fname = filename[xnn2+1:]            
        # take out .obj stufftoo
        xnn = fname.find('.obj')
        fname = fname[0:xnn] 
        tobj = bpy.data.objects[fname]
        return fname


    @property
    def name(self):
        for i in range(0,len(self.__name)):
            self.__name[i] = bpy.data.objects[self.__name[i]].name
        return self.__name

    @name.setter
    def name(self,name):
        for i in range(0,len(self.__name)):
            bpy.context.scene.objects.active = bpy.data.objects[self.__name[i]]
            bpy.data.objects[self.__name[i]].name = name[i] # object name
        self.__name = name

    @property
    def ifile(self):
        return self.__ifile

    @ifile.setter
    def ifile(self,ifile):
        from science import delete_object
        if ((ifile > 0) and (ifile is not self.__ifile)):
            # delete old files and/or surfaces
            if self.filetype[0] is not 'sphtxt':
                for n in self.name:
                    delete_object(n)
                    # now, upload a new object in your list
                    if self.filetype[0] is 'obj':
                        self.__name = self.load_obj(self.filelist[0][ifile])
                        if isinstance(self.__name, str):
                            self.__name = [self.__name]
                        self.name = self.__name
                    elif self.filetype[0] is 'ytamr':
                        #newname = [self.name + '_0']
                        ds = self.import_ytsurface(self.filelist[0][ifile],  self.isosurface_value, 
                                                   self.surf_type, self.radius, 
                                                   self.radius_units, self.surface_field, 
                                                   self.name, self.transparency, self.dist_fac, 
                                                   self.color_field, self.emit_field, self.color_map, 
                                                   self.color_log, self.emit_log, self.plot_index, 
                                                   self.color_field_max, self.color_field_min, 
                                                   self.emit_field_max, self.emit_field_min, 
                                                   self.emissivity_units)
                        self.ds = ds
                    elif self.filetype[0] is 'ytsph':
                        name, ds = self.import_ytsph(self.filelist[0][ifile], self.halo_sizes, 
                                                     self.color_field, self.color_map, 
                                                     self.color_log, self.n_ref, self.scale, 
                                                     self.coords_tag)
                        self.__name = name
                        self.name = name
                        self.ds = ds
                    else:
                        print("NOT SUPPORTED!!!!")
            else:
                # first delete everything
                #scale = self.__scale
                for n in self.name:
                    delete_object(n)
                self.__name = self.new_particle_set(self.filelist[0][ifile], self.particle_num, 
                                                    self.particle_colors, self.__halo_sizes, self.__scale)  
                self.scale = self.__scale
                self.particle_hide = self.__particle_hide
                self.name = self.__name
        self.__ifile = ifile

    @property
    def scale(self):
        for i in range(0,len(self.name)):
            self.__scale[i] = bpy.data.objects[self.name[i]].scale
        return self.__scale

    @scale.setter
    def scale(self,scale):
        if isinstance(scale,tuple):
            scale = [scale]
        for i in range(0,len(self.name)):
            bpy.context.scene.objects.active = bpy.data.objects[self.name[i]]
            bpy.data.objects[self.name[i]].scale = scale[i]
        self.__scale = scale

    @property
    def particle_hide(self):
        if self.filetype[0].find('sph') is not -1:
            for i in range(0,len(self.name)):
                self.__particle_hide[i] = bpy.data.objects[self.name[i]].hide_render
        return self.__scale

    @particle_hide.setter
    def particle_hide(self,particle_hide):
        if self.filetype[0].find('sph') is not -1:
            for i in range(0,len(self.name)):
                bpy.context.scene.objects.active = bpy.data.objects[self.name[i]]
                bpy.data.objects[self.name[i]].hide_render = particle_hide[i]
                bpy.data.objects[self.name[i]].hide = particle_hide[i]
        self.__particle_hide = particle_hide

    @property
    def halo_sizes(self):
        if self.filetype[0].find('sph') is not -1: # only have halos for sph data
            for i in range(0,len(self.name)):
                self.__halo_sizes[i] = bpy.data.materials[self.name[i]].halo.size
        return self.__halo_sizes

    @halo_sizes.setter
    def halo_sizes(self,halo_sizes):
        if len(halo_sizes) is not len(self.__halo_sizes):
            print("hey man, I need " + str(len(self.__halo_sizes)) + " numbers!")
            halo_sizes = self.__halo_sizes
        else:
            if isinstance(halo_sizes, tuple):
                hs = []
                for h in halo_sizes:
                    hs.append(h)
                halo_sizes = hs
            if self.filetype[0].find('sph') is not -1: # only have halos for sph data
                for i in range(0,len(self.name)):
                    bpy.context.scene.objects.active = bpy.data.objects[self.name[i]]
                    bpy.data.materials[self.name[i]].halo.size = halo_sizes[i]
        self.__halo_sizes = halo_sizes


    @property
    def location(self):
        for i in range(0,len(self.name)):
            self.__location[i] = bpy.data.objects[self.name[i]].location
        return self.__location

#    @location.setter
#    def location(self, pos, location):
#        if len(pos) is 0:
#            for i in range(0,len(self.name)):
#                bpy.context.scene.objects.active = bpy.data.objects[self.name[i]]
#                bpy.data.objects[self.name[i]].location = location[i]
#            self.__location = location
#        else:
#            bpy.context.scene.objects.active = bpy.data.objects[self.name[pos]]
#            bpy.data.objects[self.name[pos]].location = location
#            self.__location[pos] = location

    @location.setter
    def location(self, location):
        for i in range(0,len(self.name)):
            bpy.context.scene.objects.active = bpy.data.objects[self.name[i]]
            bpy.data.objects[self.name[i]].location = location[i]
        self.__location = location


    def export_obj(self, filename):
        # use the util one
        from scienceutils import export_obj
        export_obj(filename, self.name[0])


    def new_particle_set(self, dfile, particle_num=1, colors=None, halo_sizes=None, scale = (1,1,1)):
    # initial file, number of particles, colors array
        import bmesh
        deselect_all()
        location = (0,0,0)
        # strip "/" and ".txt"s for names
        xnn2 = 0
        xnn = 0
        # you only want the file name, not all the directory stuff
        while xnn != -1:
            xnn2 = xnn
            xnn = dfile.find('/',xnn+1)
        fname = dfile[xnn2+1:]            
        # take out .txt stufftoo
        xnn = fname.find('.txt')
        fname = fname[0:xnn] 
        pname = []
        # create different particle types
        for p in range(0,particle_num):
            #print(p)
            num = "%02d" % (p)
            pname.append(fname+'particle'+num)
            me = bpy.data.meshes.new(fname+'particleMesh'+num)
            ob = bpy.data.objects.new(fname+'particle'+num,me)
            ob.location = (0,0,0)
            bpy.context.scene.objects.link(ob)    # Link object to scene
            coords = [(0,0,0)]
            me.from_pydata(coords,[],[])
            ob.location = (0,0,0)
            if colors is not None:
                color = colors[:][p]
                halo_size = halo_sizes[p]
            else:
                print('no color')
                color = (1,1,1)
                halo_size = 0.1
            mat = makeMaterial(fname+'particle'+num, color, (1,1,1), 1.0, 1.0, mat_type='HALO', halo_size=halo_size)
            setMaterial(ob,mat) # sets everything to material 0 by default
            ob = bpy.data.objects[fname+'particle'+num] # select right object
            ob.select = True
            bpy.context.scene.objects.active=ob
        # now, add in different data types
        p2 = -1 # for toggling to object mode when using a different particle
        f = open(dfile,'r') # open file for reading
        for i, line in enumerate(f):
            pint = int(p)
            #p = data["p"][i] # get type - o - particle
            linep = line.strip() # get rid of end line
            data = linep.split(",") # split into data based on ","'s
            p = float(data[4])
            if p != p2:
                bpy.ops.object.mode_set(mode='OBJECT')  # toggle to edit mode
            num = "%02d" % (p) 
            deselect_all()
            ob = bpy.data.objects[fname+'particle'+num] # select right object
            ob.select = True
            bpy.context.scene.objects.active=ob
            bpy.ops.object.mode_set(mode='EDIT')  # toggle to edit mode
            bm = bmesh.from_edit_mesh(ob.data)
            if p !=p2: # if first instance - move origional
                if hasattr(bm.verts, "ensure_lookup_table"): # to make it work with 2.73
                    bm.verts.ensure_lookup_table()
                bm.verts[0].co = (float(data[1])*scale[pint][0],float(data[2])*scale[pint][1],float(data[3])*scale[pint][2])
                p2 = p
            else:
                bm.verts.new((float(data[1])*scale[pint][0],float(data[2])*scale[pint][1],float(data[3])*scale[pint][2]))
        f.close() # close up the file
        # update all meshes
        for p in range(0,particle_num):
            deselect_all()
            bpy.ops.object.mode_set(mode='OBJECT')  # toggle to edit mode
            num = "%02d" % (p) 
            ob = bpy.data.objects[fname+'particle'+num] # select right object
            ob.select = True
            bpy.context.scene.objects.active=ob
            bpy.ops.object.mode_set(mode='EDIT')  # toggle to edit mode
            bmesh.update_edit_mesh(ob.data)
        bpy.ops.object.mode_set(mode='OBJECT')  # toggle to edit mode
        return pname


    # internal importer from yt
    def _import_ytsurface(self, vertices, colors, alpha, emisses, colorindex, 
                          meshname='Steve', meshlocation=(0,0,0), plot_index=0):
        import numpy as np
        deselect_all()
        nftype = [("f1","int"),("f2","int"),("f3","int"),("f4","int")]
        newfaces = np.empty(int(len(vertices)/3.), dtype=nftype) # store sets of face colors
        cc = 0
        for i in range(0,int(len(vertices)/3.)):
            newfaces[i] = (cc, cc+1, cc+2, cc) # repeat last for triangles
            cc = cc+3
        me = bpy.data.meshes.new(meshname+"Mesh")
        ob = bpy.data.objects.new(meshname,me)
        ob.location = meshlocation   # position object at 3d-cursor
        bpy.context.scene.objects.link(ob)                # Link object to scene
        # Fill the mesh with verts, edges, faces 
        me.from_pydata(vertices.tolist(),[],newfaces.tolist())
        me.update(calc_edges=True)    # Update mesh with new data
        # materials woop woop!
        obj = bpy.data.objects[meshname]
        bpy.context.scene.objects.active=obj
        for i in range(0,len(colors[0])):
            mat = makeMaterial(meshname+str(i)+'_'+str(plot_index), 
                               (colors[0][i],colors[1][i],colors[2][i]), 
                               (1,1,1), alpha, emisses[i])
            setMaterial(obj,mat)
        # now, do individual mat indecies
        for i in range(0,len(newfaces)):
            me.polygons[i].material_index = colorindex[i]


    def import_ytsurface(self, filename, isosurface_value = 1e-27, surf_type='sphere', radius = 10.0, 
                         radius_units = "mpc", surface_field="density",  
                         meshname = 'Allen',  
                         transparency = 1.0, dist_fac = None,
                         color_field = None, emit_field = None, color_map = "algae", 
                         color_log = True, emit_log = True, plot_index = None, 
                         color_field_max = None, color_field_min = None, 
                         emit_field_max = None, emit_field_min = None, emissivity_units = None):
        from yt import load
        import numpy as np
        pf = load(filename)

        # for only 1 value
        if isinstance(isosurface_value, float):
            isosurface_value = [isosurface_value]
            transparency = [transparency]
            plot_index = 0
        if surf_type is 'sphere':
            presurf = pf.sphere("max", (radius, radius_units)) 
        else:
            print("Waaaahhoooo, don't know anything else!  Am doing spherical surface for the time being")
            presurf = pf.sphere("max", (radius, radius_units)) 
        if emit_field is not None:
            if emissivity_units is None:
                emissivity_units = "g**2*K**0.5/cm**6"
            pf.add_field(("gas","emissivity"), units=emissivity_units, function = emit_field)
            emit_field = 'emissivity'
        for i in range(0,len(isosurface_value)):
            trans = transparency[i]
            surf = pf.surface(presurf, surface_field, isosurface_value[i])
            m_name = meshname + '_' + str(i)
            vertices, colors, alpha, emisses, colorindex = surf.export_blender(transparency = trans, dist_fac = dist_fac,
                                                                               color_field = color_field, emit_field = emit_field, color_map = color_map, 
                                                                               color_log = color_log, emit_log = emit_log, plot_index = i, 
                                                                               color_field_max = color_field_max, color_field_min = color_field_min, 
                                                                               emit_field_max = emit_field_max, emit_field_min = emit_field_min)
            self._import_ytsurface(vertices, colors, alpha, emisses, colorindex, m_name, (0,0,0), i)
        return pf



    # import sph data via yt
    def import_ytsph(self, filename, halo_sizes, color_field, color_map, color_log, n_ref, scale, coords_tag):

        from yt import load
        import numpy as np
        import bmesh

        ds = load(filename)

        # strip "/" for naming
        xnn2 = 0
        xnn = 0
        # you only want the file name, not all the directory stuff
        while xnn != -1:
            xnn2 = xnn
            xnn = filename.find('/',xnn+1)

        fname_out = filename[xnn2+1:]            

        # get coords and color data
        dd = ds.all_data()
        xcoord = dd[coords_tag][:,0].v
        ycoord = dd[coords_tag][:,1].v
        zcoord = dd[coords_tag][:,2].v
        cs = dd[color_field]

        # map colorcode to 256 material colors
        if color_log: cs = np.log10(cs)

        mi, ma = cs.min(), cs.max()
        cs = (cs - mi) / (ma - mi)

        from yt.visualization._colormap_data import color_map_luts # import colors 
        # the rgb colors
        colors = color_map_luts[color_map]

        x = np.mgrid[0.0:1.0:colors[0].shape[0]*1j]
        # how the values map to the colors
        color_index = (np.interp(cs,x,x)*(colors[0].shape[0]-1)).astype("uint8")
        color_index_list = color_index.tolist()

        name = []
        hused = []
        # create scale if necessary
        if len(scale[:][:]) < color_index.max(): # generate larger scale
            sc = (scale[0][0], scale[0][1], scale[0][2])
            scale = []
            for i in range(color_index.min(), color_index.max()):
                scale.append(sc)

        # also, allow for multiple halo sizes
        if len(halo_sizes) < color_index.max():
            hs = halo_sizes[0]
            halo_sizes = []
            for i in range(color_index.min(), color_index.max()):
                halo_sizes.append(hs)

        # create sph data into meshes and color them
        for ind in range(color_index.min(), color_index.max()):
            cl = [i for i,j in enumerate(color_index_list) if j == ind]
            if len(cl) > 0:
                fname = fname_out + '_' + str(ind)
                name.append('particle_'+fname)
                hused.append(halo_sizes[ind])
                me = bpy.data.meshes.new('particleMesh_'+fname)
                ob = bpy.data.objects.new('particle_'+fname,me)
                ob.location = (0,0,0)
                bpy.context.scene.objects.link(ob)    # Link object to scene
                coords = [(0,0,0)]
                me.from_pydata(coords,[],[])
                ob.location = (0,0,0)
                ob = bpy.data.objects['particle_'+fname] # select right object
                deselect_all()
                ob.select = True
                bpy.context.scene.objects.active=ob
                mat = makeMaterial('particle_'+fname, 
                                           (colors[0][ind],colors[1][ind],colors[2][ind]), 
                                           (1,1,1), 1.0, 1.0, mat_type = 'HALO', halo_size=halo_sizes[ind])
                setMaterial(ob,mat)
                # add in verts
                bpy.ops.object.mode_set(mode='EDIT')  # toggle to edit mode
                bm = bmesh.from_edit_mesh(ob.data)
                # now, find all verts with this color index (ind)
                # move original vertex to actual location
                if hasattr(bm.verts, "ensure_lookup_table"): # to make it work with 2.73
                    bm.verts.ensure_lookup_table()
                bm.verts[0].co = (xcoord[cl[0]]*scale[ind][0],ycoord[cl[0]]*scale[ind][1],zcoord[cl[0]]*scale[ind][2])
                for i in range(1,len(xcoord[cl])):
                    bm.verts.new((xcoord[cl[i]]*scale[ind][0],ycoord[cl[i]]*scale[ind][1],zcoord[cl[i]]*scale[ind][2]))
                bmesh.update_edit_mesh(ob.data)
                bpy.ops.object.mode_set(mode='OBJECT')
        return name, ds, hused
