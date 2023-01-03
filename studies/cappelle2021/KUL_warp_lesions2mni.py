#!/usr/bin/env python

from dis import dis
import glob 
import os


bidsdir = './T1T2FLAIRMTR_ratio'
outdir = './samseg_long'
"""
for root, dirs, files in os.walk(bidsdir):

    #print(dirs)
    #print(files)

    for dir in dirs:

        if 'ses-' in dir: 
            if not 'fs' in dir: 
                #print(dirs)

                searchdir = os.path.join(root, dir)
                #print(searchdir)


                ImType = "T1w"
                
                searchIms = searchdir + '/*space-MNI_' + ImType + '.nii.gz'
                #print(searchIms)
                
                # find all Im
                Ims = glob.glob(searchIms, recursive=False)
                #print(Ims)
                    

                for Im in Ims:
                    #print(Im)
                    # find the warp file
                    # warp lesions to MNI
                    ses = os.path.split(os.path.splitext(os.path.splitext(Im)[0])[0])[0]
                    base_name1 = os.path.split(os.path.splitext(os.path.splitext(Im)[0])[0])[1]
                    base_name = base_name1.split('_space-MNI_T1w')[0]
                    #print(ses)
                    #print(base_name)
                    transform1 = os.path.join(ses,'warp2mni',base_name) + '_T1w2MNI_1Warp.nii.gz'
                    transform2 = '[' + os.path.join(ses,'warp2mni',base_name) + '_T1w2MNI_0GenericAffine.mat' + ',0]'
                    #print(transform1)
                    #print(transform2)
                    input = os.path.join(ses,'rois',base_name) + '_MSLesion.nii.gz'
                    output = os.path.join(ses,base_name) + '_space-MNI_MSLesion.nii.gz'
                    #print(input)
                    #print(output)
                    
                    #antsApplyTransforms -d 3 \
                    #    --verbose $ants_verbose \
                    #    -i $input \
                    #    -o $output \
                    #    -r $reference \
                    #    -t $transform1 -t $transform2 \
                    #    -n $interpolation_type
                    cmd = 'antsApplyTransforms -d 3  -i ' + input + \
                        ' -o ' + output + \
                        ' -r ' + Im + \
                        ' -t ' + transform1 + ' ' + transform2 + \
                        ' -n NearestNeighbor'

                    #print(cmd)
                    
                    if not os.path.exists(output):
                        print('Warping image ' + input)
                        out = os.popen(cmd).read().strip()
                        print(out)
                    else:
                        print('Already warped ' + input)
"""
# Compute the sum of lesions of the first timepoint of all subjects
old_root=[]
list = []
for root, dirs, files in os.walk(bidsdir):    
    for dir in dirs:
        if 'ses-' in dir and not 'fs' in root: 
            if not old_root == root:

                if len(dirs) > 1: 
                    im = os.path.join(root,dirs[0],'*MSLesion.nii.gz')
                    list.append(im)
            #else:
            #    print('gelijk')
            old_root=root
#print(' '.join(list))
print(len(list))
cmd = 'mrmath ' + ' '.join(list) + ' sum ' + 'tp_first_MSLesion.nii.gz -force'
out = os.popen(cmd).read().strip()
print(out)

# Compute the sum of lesions of the last timepoint of all subjects
# Note if a subject has only 1 session, this is also included (as in first)       
old_root=[]
list = []
for root, dirs, files in os.walk(bidsdir):    
    for dir in dirs:
        if 'ses-' in dir and not 'fs' in root: 
            if not old_root == root:
                if len(dirs) > 1: 
                    im = os.path.join(root,dirs[-1],'*MSLesion.nii.gz')
                    list.append(im)
            #else:
            #    print('gelijk')
            old_root=root
#print(' '.join(list))
print(len(list))
cmd = 'mrmath ' + ' '.join(list) + ' sum ' + 'tp_last_MSLesion.nii.gz -force'
out = os.popen(cmd).read().strip()
print(out)             
