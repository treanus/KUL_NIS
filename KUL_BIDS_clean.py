#!/usr/bin/env python

import glob 
import os
import numpy as np
import json

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
                elif nIms > 1:
                    if ImType == "T2w":
                        # we need to keep the transverse
                        print('There are ' + str(nIms) + ' ' + ImType + ' images')
                        print ('Notably:')
                        i=0
                        orientation = np.zeros((nIms,3))
                        for Im in Ims:
                            print(Im)
                            ImJson = os.path.splitext(os.path.splitext(Im)[0])[0] + '.json'
                            #print(ImJson)
                            with open(ImJson, 'r') as myfile:
                                data=myfile.read()
                            obj = json.loads(data)
                            orientationfull = obj['ImageOrientationPatientDICOM']
                            #print(orientationfull)
                            orientation[i] = orientationfull[0:3]
                            i = i + 1 
                        #print(orientation)
                        ori=np.argmax(orientation, axis=1)
                        #print(ori)
                        keep = np.argmin(ori)
                        print('Keeping ' + Ims[keep])
                        
                    elif ImType == "T1w":
                        # we need to keep the youngest 3D
                        print('There are ' + str(nIms) + ' ' + ImType + ' images')
                        print ('Notably:')
                        i=0
                        seriesnum = np.zeros((nIms,1))
                        acq = []
                        spacing = np.zeros( (nIms,3))
                        for Im in Ims:
                            print(Im)
                            # read seriesnum
                            ImJson = os.path.splitext(os.path.splitext(Im)[0])[0] + '.json'
                            #print(ImJson)
                            with open(ImJson, 'r') as myfile:
                                data=myfile.read()
                            obj = json.loads(data)
                            seriesnum[i] = obj['SeriesNumber']
                            # read acquisitiontype
                            acq.append(obj['MRAcquisitionType'])
                            i = i + 1 
                        #print(seriesnum)
                        #print(acq)
                        if '3D' in acq:
                            keep_3D = list(filter(lambda i: acq[i]=="3D", range(len(acq))))
                            filtered_seriesnum = seriesnum[keep_3D]
                        else:
                            filtered_seriesnum = seriesnum
                        #print(filtered_seriesnum)
                        msn = np.min(filtered_seriesnum) #minimal seriesnumbr (youngest)
                        #print(msn)
                        keep = np.where(seriesnum == msn)
                        keep = np.asscalar(keep[0])
                        #print(keep)
                        print('Keeping ' + Ims[keep])
                        
                    elif ImType == "FLAIR":
                        # we need to keep the highest resolution
                        keep=0
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
                        voxelvolume = np.prod(spacing,axis=1)
                        print(voxelvolume)
                        #maxvoxelsize = np.max(spacing, axis=1)
                        keep=np.argmin(voxelvolume)
                        print('Keeping ' + Ims[keep])

                    # Now do the change
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
