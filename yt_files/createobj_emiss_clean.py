import yt
import numpy as np
    
# different data - FLASH
filename = '~/data/GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100'
sphere_rad = 15.0 # in kpc

#Enzo
#filename = '~/data/IsolatedGalaxy/galaxy0030/galaxy0030'
#sphere_rad = 200.0 # in kpc


outfile = 'galsurfaces'

rho = 1e-27 # for each surface
trans = 1.0 # for transparency of each surface

color_field = 'temperature' # color your surface by this

pf = yt.load(filename)


# emissivity of the material
# this needs to be a combination of the color_field and surface field
#def _Emissivity(field, data):
#    return (data['gas','density']*data['density']*np.sqrt(data['gas','temperature']))

#pf.add_field(("gas","emissivity"), units="g**2*K**0.5/cm**6", function=_Emissivity)


dd = pf.h.sphere("max", (sphere_rad, "kpc"))

surf = pf.h.surface(dd, 'density', rho)
surf.export_obj(outfile, transparency = trans, 
                color_field=color_field, emit_field = 'emissivity')
        
