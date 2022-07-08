#!/usr/bin/env python

import argparse
import glob 
import os
import numpy as np
import nibabel as nib
from scipy import ndimage

parser = argparse.ArgumentParser(description="Determine the postion of the DRT on AC-PC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-v", "--verbose", action="store_true", help="increase verbosity")
parser.add_argument("-n", "--ncpu", help="number of threads to use for samseg")
parser.add_argument("dest", help="Destination location")
args = parser.parse_args()
config = vars(args)
#print(config)

if args.ncpu is None:
    ncpu = 15
else:
    ncpu = args.ncpu

bidsdir = './fmriprep'
outdir = args.dest
#print(outdir)

results_csv = os.path.join(outdir) + 'Results_DRT.csv'
cmd = 'echo "base_name, type, ses, side, count, CMx, CMy, CMz" > ' + results_csv
print(cmd)
out = os.popen(cmd).read().strip()
print(out)

for root, dirs, files in os.walk(bidsdir):
    for dir in dirs:
        if 'sub-' in dir:
            searchdir = os.path.join(root, dir)
            #print(searchdir)
            #print(dir)
            os.makedirs(os.path.join(outdir,dir), exist_ok=True)

            for ImType in ["space-MNI152NLin2009cAsym_desc-preproc_T1w"]:

                searchIms = searchdir + '/anat/*' + ImType + '.nii.gz'
                print(searchIms)
                

                # find all Im
                Ims = glob.glob(searchIms)
                print(str(Ims))

                for Im in Ims: 
                    # empty the MNI image, except slice 95 (AC/PC)
                    dir_name, base_name = os.path.split(os.path.splitext(os.path.splitext(Im)[0])[0])
                    output = os.path.join(outdir, dir, base_name) + '_acpc_plane.nii.gz'
                    print(output)
                    cmd = 'mrgrid ' + Im + ' crop -axis 2 95,144 - | mrgrid - pad -axis 2 95,144 -force ' + output
                    print(cmd)
                    out = os.popen(cmd).read().strip()
                    print(out)
                    
                    cmd = 'cp ' + Im + ' ' + os.path.join(outdir,dir)
                    out = os.popen(cmd).read().strip()
                    print(out)

                    # warp that back to subject space
                    input = output
                    base_name = base_name.split('_space')[0]
                    print(base_name)
                    plane = os.path.join(outdir, dir, base_name) + '_T1w_acpc_plane.nii.gz'
                    transform = os.path.join(dir_name, base_name) + '_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5'
                    reference = os.path.join(dir_name, base_name) + '_desc-preproc_T1w.nii.gz'
                    cmd = 'antsApplyTransforms -d 3 --float 1 --verbose 1' + \
                        ' -i ' + input + \
                        ' -o ' + plane + \
                        ' -r ' + reference + \
                        ' -t ' + transform + \
                        ' -n Linear'
                    print(cmd)
                    out = os.popen(cmd).read().strip()
                    print(out)
                    cmd = 'cp ' + reference + ' ' + os.path.join(outdir,dir)
                    out = os.popen(cmd).read().strip()
                    print(out)

                    if 0 : 
                        # intersect the DRT with the plane in subject space
                        drt = os.path.join('.','BIDS','derivatives','KUL_compute',base_name,'ses-T0','FWT', base_name + '_TCKs_output','DRT_LT_output', 'DRT_LT_fin_BT_iFOD2.tck')
                        print(drt)
                        cmd = 'cp ' + drt + ' ' + os.path.join(outdir,dir)
                        print(cmd)
                        out = os.popen(cmd).read().strip()
                        print(out)

                        # smooth the tract
                        drt_smooth = os.path.join(outdir, dir, base_name) + '_DRT_LT_smooth.tck'
                        cmd = 'scil_smooth_streamlines.py -f --gaussian 25 --reference ' + plane + ' ' + drt + ' ' + drt_smooth
                        out = os.popen(cmd).read().strip()
                        print(out)
                        out = os.popen(cmd).read().strip()
                        print(out)

                        # make a tckmap
                        drt_map = os.path.join(outdir, dir, base_name) + '_DRT_LT_smooth_map.nii.gz'
                        cmd = 'tckmap -force -contrast tdi -template ' + plane + ' ' + drt_smooth + ' ' + drt_map
                        print(cmd)
                        out = os.popen(cmd).read().strip()
                        print(out)
                    
                    i = 1

                    for side in ['LT', 'RT']:
                        for ses in ['T0','T1','T2']:

                            drt = os.path.join('.','BIDS','derivatives','KUL_compute',base_name,'ses-' + ses,'FWT', base_name + \
                                '_TCKs_output','DRT_' + side + '_output', 'DRT_' + side + '_fin_map_BT_iFOD2.nii.gz')
                            if os.path.exists(drt):
                                
                                print(drt)
                                cmd = 'cp ' + drt + ' ' + os.path.join(outdir,dir)
                                print(cmd)
                                out = os.popen(cmd).read().strip()
                                print(out)
                                
                                # find the number of streamlines
                                tck = os.path.join('.','BIDS','derivatives','KUL_compute',base_name,'ses-' + ses,'FWT', base_name + \
                                '_TCKs_output','DRT_' + side + '_output', 'DRT_' + side + '_fin_BT_iFOD2.tck')
                                cmd = 'tckstats -output count ' + tck 
                                print(cmd)
                                out = os.popen(cmd).read().strip()
                                count = out.splitlines()[0]
                                print(count)

                                # regrid to HR T1W
                                regrid = os.path.join(outdir, dir, base_name) + '_DRT_' + side + '_ses-' + ses + '_map_regrid.nii.gz'
                                cmd = 'mrgrid -force ' + drt + ' regrid -template ' + plane + ' ' + regrid 
                                print(cmd)
                                out = os.popen(cmd).read().strip()
                                print(out)

                                # make an intersection image
                                intersect = os.path.join(outdir, dir, base_name) + '_DRT_' + side + \
                                    '_ses-' + ses + '_map_intersect.nii.gz'
                                cmd = 'mrcalc -force ' + plane + ' 1 -gt ' + regrid + ' -mul ' + intersect
                                print(cmd)
                                out = os.popen(cmd).read().strip()
                                print(out)

                                # find the center of mass and write as an image
                                img = nib.load(intersect)
                                img_data = img.get_data()
                                CM = ndimage.measurements.center_of_mass(img_data)
                                print(CM)
                                #print(round(CM[0]))
                                cmd = 'echo "' + base_name + ', CM,' + ses + ',' + side + ',' + count + ',' + \
                                    str(CM[0]) + ',' +  str(CM[1]) + ',' + str(CM[2]) + '" >> ' + results_csv
                                print(cmd)
                                out = os.popen(cmd).read().strip()
                                print(out)
                                
                                #print(img.header.get_data_shape())
                                CM_image = os.path.join(outdir, dir, base_name) + '_DRT_' + side + \
                                    '_ses-' + ses + '_cm.nii.gz'
                                new_img_data = np.zeros(img.header.get_data_shape())
                                new_img_data[round(CM[0]), round(CM[1]), round(CM[2])] = i

                                new_img = nib.Nifti1Image(new_img_data, img.affine, img.header)
                                nib.save(new_img, CM_image)
                                

                                # find the voxel with most streamlines
                                mp = ndimage.measurements.maximum_position(img_data)
                                print(mp)
                                cmd = 'echo "' + base_name + ', mp,' + ses + ',' + side + ',' + count + ',' + \
                                    str(mp[0]) + ',' +  str(mp[1]) + ',' + str(mp[2]) + '" >> ' + results_csv
                                print(cmd)
                                out = os.popen(cmd).read().strip()
                                print(out)

                                mp_image = os.path.join(outdir, dir, base_name) + '_DRT_' + side + \
                                    '_ses-' + ses + '_mp.nii.gz'
                                new_img_data = np.zeros(img.header.get_data_shape())
                                new_img_data[round(mp[0]), round(mp[1]), round(mp[2])] = i

                                new_img = nib.Nifti1Image(new_img_data, img.affine, img.header)
                                nib.save(new_img, mp_image)
                                i = i + 1

                                # make an outline
                                outline = os.path.join(outdir, dir, base_name) + '_DRT_' + side + \
                                    '_ses-' + ses + '_map_intersect_outline.nii.gz'
                                cmd = 'maskfilter ' + intersect + ' dilate - | mrcalc -force - ' + intersect + ' -sub ' + outline
                                print(cmd)
                                out = os.popen(cmd).read().strip()
                                print(out)
                            
                            else: 
                                print('No DRT found!')
                                cmd = 'echo "' + base_name + ',' + ses + ',' + side + ',' + 'NaN,NaN, NaN, NaN" >> ' + results_csv
                                print(cmd)
                                out = os.popen(cmd).read().strip()
                                print(out)
                
