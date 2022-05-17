#!/usr/bin/env python

import glob 
import os
import numpy as np

bidsdir = './BIDS'

for subdir, dirs, files in os.walk(bidsdir):
    for dir in dirs:
        if 'anat' in dir:
            searchdir = os.path.join(subdir, dir)
            #print(searchdir)

            for ImType in ["T1w", "T2w", "FLAIR"]:

                searchIms = searchdir + '/*' + ImType + '.nii.gz'
                #print(searchIms)

                # find all Im
                Ims = glob.glob(searchIms)
                nIms = len(Ims)

                if nIms == 0:
                    print('No ' + ImType + ' images, doing nothing')
                elif nIms == 1:
                    print('There is only one ' + ImType + ', keeping this one')
                else:
                    print('There are ' + str(nIms) + ' ' + ImType + ' images')
                    spacing = np.zeros( (nIms,3))
                    print ('Notably:')
                    i=0
                    for Im in Ims:
                        print(Im)
                        cmd = 'mrinfo ' + Im + ' -spacing'
                        spacing[i] = os.popen(cmd).read().strip().split()
                        i = i + 1    
                    print(spacing)
                    maxvoxelsize = np.max(spacing, axis=1)
                    keep=np.argmin(maxvoxelsize)
                    print('Keeping ' + Ims[keep])
                    i = 0
                    for Im in Ims:
                        if i == keep:
                            p1 = Im.split('_run')[0]
                            #print(p1)
                            cmd1 = 'mv ' + Im + ' ' + p1 + '_' + ImType + '.nii.gz'
                            p3 = Im.split('.nii.gz')[0] + '.json'
                            cmd2 = 'mv ' + p3 + ' ' + p1 + '_' + ImType + '.json'
                        else:
                            cmd1 = 'rm -f ' + Im
                            p1 = cmd1.split('.nii.gz')[0]
                            cmd2 = p1 + '.json'
                        print(cmd1)
                        print(cmd2)
                        
                        out = os.popen(cmd1).read().strip()
                        print(out)
                        out = os.popen(cmd2).read().strip()
                        print(out)
                        i = i + 1
