    #Copyright 2009-2016 Seyed Hessam Moosavi Mehr, Juergen Probst
    #This program is free software; you can redistribute it and/or modify
    #it under the terms of the GNU General Public License as published by
    #the Free Software Foundation; either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.

    #You should have received a copy of the GNU General Public License
    #along with this program. If not, see <http://www.gnu.org/licenses/>.

speed_of_light = 299792458 # m/s

# for 1500nm: ? - (from first submitted data.py)
#dielectrics = {'si':11.8,'diamond':5.7,'pe':2.26,\
#'sio2':4.5,'paper':3.5,'teflon':2.1,'pvc':4.5}

# for 650nm
wlen=650
refr_index = {
    'SiN'     : 2.0085,
    '3C-SiC'  : 2.55378 + 3.417e4 / wlen**2, # Schaffer 1969 (http://www.ioffe.ru/SVA/NSM/Semicond/SiC/optic.html)
    # ordinary direction: polarisation perpendicular to the optic axis of material:
    '4H-SiC-o': 2.5610 + 3.40e4 / wlen**2, 
    # extraordinary direction: polarisation parallel to the optic axis of material:
    '4H-SiC-e': 2.6041 + 3.75e4 / wlen**2,
    '6H-SiC-o': 2.5531 + 3.34e4 / wlen**2, 
    '6H-SiC-e': 2.5852 + 3.68e4 / wlen**2,
    # source 4H & 6H: Schaffer1971 https://www.osapublishing.org/ao/abstract.cfm?uri=ao-10-5-1034
    'SiO2'    : 1.5,
    'PMMA950' : 1.4957, # Microchem 950 data sheet
}
# optic axis of material in z direction:
refr_index['4H-SiC-anisotropic_c_in_z'] = (
    refr_index['4H-SiC-o'], 
    refr_index['4H-SiC-o'], 
    refr_index['4H-SiC-e'])
refr_index['6H-SiC-anisotropic_c_in_z'] = (
    refr_index['6H-SiC-o'], 
    refr_index['6H-SiC-o'], 
    refr_index['6H-SiC-e'])

# short material names for use in graphics titles:    
material_names = dict([(k, k) for k in refr_index.keys()])
material_names['4H-SiC-anisotropic_c_in_z'] = '4H-SiC-z'
material_names['6H-SiC-anisotropic_c_in_z'] = '6H-SiC-z'
    
    
#print 'refractive index:', refr_index

dielectrics = dict(refr_index)
for key, val in dielectrics.items():
    if isinstance(val, (tuple, list)):
        dielectrics[key] = tuple(v**2 for v in val)
    else:
        dielectrics[key]=val**2

#print 'epsilon:', dielectrics
    



