import science

#tipsy data
filename = '/Users/jillnaiman/data/TipsyGalaxy/galaxy.00300'
#color_field = ('Gas', 'Temperature') 
#color_log = True
#color_map = 'Rainbow'
#scale = [(1.0, 1.0, 1.0)]


# Gadget data
filename = '/Users/jillnaiman/data/GadgetDiskGalaxy/snapshot_200.hdf5'
scale = [(0.1, 0.1, 0.1)]

# these two things play off eachother!
halo_size = 0.108 # need to play with this
set_cam = (0,0,70)


cam = science.Camera()
cam.location = set_cam
cam.clip_begin = 0.0001

lighting = science.Lighting('EMISSION')

# initialize render
render_directory = '/Users/jillnaiman/blenderRenders/'
render_name = 'gadgetRender_'
render = science.Render(render_directory, render_name)

#tipsy
#myobject = science.Load(filename, scale=scale, halo_sizes = halo_size, 
#                        color_map = color_map, 
#                        color_log = color_log, n_ref=8)

#gadget
myobject = science.Load(filename, scale=scale, halo_sizes = halo_size)

#render.render()
