#!/usr/bin/env python3

# https://neurostars.org/t/extract-voxel-coordinates/7282
# http://blog.chrisgorgolewski.org/2014/12/how-to-convert-between-voxel-and-mm.html
# https://stackoverflow.com/questions/6967463/iterating-over-a-numpy-array
# https://stackabuse.com/calculating-euclidean-distance-with-numpy/
# https://numpy.org/doc/stable/reference/arrays.nditer.html
# https://thispointer.com/numpy-amin-find-minimum-value-in-numpy-array-and-its-index/

# Author: Ahmed Radwan, UZ Leuven/KU Leuven - ahmed.radwan@kuleuven.be, radwanphd@gmail.com
# This python script was developed in python3.8
# v: 0.6 - 09122021
# This workflow does the following:
# 1- Recieve input (described in def main lines 23-47) from user 
# 2- Find the external edge of the 2 input masks
# 3- Find centers of gravity from numpy arrays of both input masks 
# 4- Calculate Euclidean distances between every voxel of A to every voxel of B
# 5- Calculate distances between every voxel of A to cogB and vice versa
# 6- Calculate distances between both masks' COGs
# 7- Find index of voxels giving min distances
# 8- Print to CLI and save to output text files and nifti files
# 9- This is supplemented by overlap COG, respective distance calculations and overlap count and volume ratios if masks are initially overlapping

## To do:
# 1- What if more than one voxel have shortest distance?

import os, sys, getopt
import nibabel as nib
import numpy as np
from scipy import ndimage

# define main input function here1
def main(argv):
    ilog = ''
    inii = ''
    iname = ''
    ofolder = ''
    try:
        opts, args = getopt.getopt(argv,"ha:b:o:",["in1=","in2=","o="])
    except getopt.GetoptError:
        print ('KUL_EDs_b2masks.py -a <in1> -b <in2> -o <out>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('KUL_EDs_b2masks.py calculates Euclidean distances between two binary masks')
            print ('KUL_EDs_b2masks.py will also check for initial overlap and calculate distances to and from overlapping voxels as well')
            print ('The two input masks must be in the same space and have the same dimensions')
            print ('The first mask should be the smaller one (e.g. DES sphere, or lesion mask), and the second the larger (e.g. CST)')
            print ('KUL_EDs_between_2masks.py -a <in1> -b <in2> -o <out>')
            sys.exit()
        elif opt in ("-a", "--in1"):
            in1 = arg
        elif opt in ("-b", "--in2"):
            in2 = arg
        elif opt in ("-o", "--out"):
            out = arg
    print ('Input full path and file name for the first mask image "', in1)
    print ('Input full path and file name for the second mask image "', in2)
    print ('Prefix output name "', out)

    # for debugging
    # in1 = '/media/radwan/AR_16T/S61759_BIDS_fMRI/BIDS/derivatives/Warping_2_native/ECS/sub-PT004_ECS2nat/sub-PT004_ECS_split/Spheres_split_2_reconned.nii.gz'
    # in2 = '/media/radwan/AR_16T/S61759_BIDS_fMRI/BIDS/derivatives/Warping_2_native/ECS/sub-PT004_ECS2nat/sub-PT004_ECS_split/Spheres_split_3_reconned.nii.gz'
    # in2 = '/media/radwan/AR_16T/S61759_BIDS_fMRI/BIDS/derivatives/Warping_2_native/TCKs/sub-PT004_TCKs_warping/TCK_maps/AF_all_all_LT_fin_BT_map_inNat.nii.gz'
    # out = 'Alpha_trial'

    # now we load in the niis
    img1 = nib.load(in1)
    img2 = nib.load(in2)

    # grab their affines
    aff1 = img1.affine
    aff2 = img2.affine

    # sanity check, are the affines the same or close enough ?
    if np.allclose(aff1, aff2):
        # for each input convert fdata to 16bit uint
        # then do a 1x iterative cleaning, 1x erosion, absolute difference is the edge
        im1_data = np.uint16(img1.get_fdata())
        im2_data = np.uint16(img2.get_fdata())

        # here we start checking for initial overlaps
        # if this is found we can follow a different workflow
        # where the overlap voxels COG is calculated
        # and we continue while focusing on those instead of the whole image
        in_overlap = np.uint16(np.multiply(im1_data, im2_data))

        # if the initial overlap is zero we do the wf described above
        # if not we look at the overlapping voxels
        # mask B voxels overlapping with mask A
        # get their COG
        # calculate distance between COGA and COG of maskB voxels overlapping with maskA (and vice versa?)
        # also calculate percentage of maskB overlapping with maskA and vice versa
        if np.amax(in_overlap) != 0:
            wf = 1
            print('overlap found')
        else:
            wf = 2
            print('overlap not found')

        # clean_im1 = np.uint16(ndimage.morphology.binary_dilation((ndimage.morphology.binary_erosion \
        #     (im1_data, iterations=2)), iterations=1))
        
        eroded_im1 = np.uint16(ndimage.morphology.binary_erosion(im1_data))
        eroded_im2 = np.uint16(ndimage.morphology.binary_erosion(im2_data))

        # clean_im2 = np.uint16(ndimage.morphology.binary_dilation((ndimage.morphology.binary_erosion \
        #     (im2_data, iterations=2)), iterations=1))
        
        outline1 = np.uint16(np.absolute(np.subtract(im1_data,eroded_im1)))
        outline2 = np.uint16(np.absolute(np.subtract(im2_data,eroded_im2)))
        
        # to get indices of nonzero voxels
        img1_idx = np.where(outline1)
        img2_idx = np.where(outline2)

        # to get cogs in voxel coords
        cog1 = ndimage.measurements.center_of_mass(im1_data)
        cog2 = ndimage.measurements.center_of_mass(im2_data)
        # then convert voxel coords to mm
        cog1_xyz = nib.affines.apply_affine(aff1, cog1)
        cog2_xyz = nib.affines.apply_affine(aff2, cog2)

        # list of arrays to (voxels, 3) array
        ijk1 = np.vstack(img1_idx).T
        ijk2 = np.vstack(img2_idx).T

        # convert the voxel coordinates to mm coordinates
        xyz1 = nib.affines.apply_affine(aff1, ijk1)
        xyz2 = nib.affines.apply_affine(aff2, ijk2)

        # declare empty numpy arrays
        # to enable recovery of voxel coordinates afterwards
        results = np.zeros((ijk1.shape[0], ijk2.shape[0]), np.float32)
        cog1_ds = np.zeros((ijk2.shape[0]), np.float32)
        cog2_ds = np.zeros((ijk1.shape[0]), np.float32)

        vox_A_maps = np.zeros(im1_data.shape, np.uint16)
        vox_B_maps = np.zeros(im2_data.shape, np.uint16)
        COGA_map = np.zeros(im1_data.shape, np.uint16)
        COGB_map = np.zeros(im2_data.shape, np.uint16)

        # what is the distance between the COGs of both masks
        cogs_d = np.linalg.norm(cog1_xyz-cog2_xyz)

        # loop 1 over in1 nonzero voxel mm coordinates
        # calculate cog2 distance to every voxel in in1
        # loop 2 over in2 nonzero voxel mm coordinates
        # calculate distances between every voxel in in1 to in2
        for ii in range(0,xyz1.shape[0]):
            cog2_ds[ii] = np.linalg.norm(cog2_xyz-xyz1[ii])
            for jj in range(0,xyz2.shape[0]):
                results[ii,jj] = np.linalg.norm(xyz1[ii]-xyz2[jj])
                

        # loop 3 over nonzero voxels in in2
        # calculate distances between cog1 and every voxel in in2
        for uu in range(0,xyz2.shape[0]):
            cog1_ds[uu] = np.linalg.norm(cog1_xyz-xyz2[uu])


        # convert all to numpy arrays for safety
        # find min ds
        all_min = (np.amin(results))
        coga_2b = (np.amin(cog1_ds))
        cogb_2a = (np.amin(cog2_ds))

        # find index of min distance entry
        alidx = np.where(results == all_min)
        
        # grab the coordinates of the voxels giving shortest ds from both masks
        a_vox_mm = xyz1[alidx[0]]
        a_vox_vv = ijk1[alidx[0]]
        b_vox_mm = xyz2[alidx[1]]
        b_vox_vv = ijk2[alidx[1]]
        
        # find coordinates in vox and mm of the voxels with shortest distances
        ca2bijk = ijk2[np.where(cog1_ds == coga_2b)]
        ca2bxyz = xyz2[np.where(cog1_ds == coga_2b)]
        cb2aijk = ijk1[np.where(cog2_ds == cogb_2a)]
        cb2axyz = xyz1[np.where(cog2_ds == cogb_2a)]

        print('COG of mask A in mm:', cog1_xyz, 'in voxels: ', cog1)
        print('COG of mask B in mm:', cog2_xyz, 'in voxels: ', cog2)
        print('Minimum distance between all voxels of mask A and all voxels of mask B = ', all_min)
        print('Minimum distance between COG of mask A and all voxels of mask B = ', coga_2b)
        print('Minimum distance between COG of mask B and all voxels of mask A = ', cogb_2a)
        print('Minimum distance between COG of mask A and COG of mask B = ', (cogs_d))

        # define output dirs and files
        pwd = str(os.popen('pwd').read()).strip()

        # is the output given?
        if 'out' in locals():
            out_n = out
        else:
            out_n = 'KUL_EDs'


        # make output dir
        os.system('mkdir -p ' + pwd + '/' + out_n + '_output')

        # get basename of input files
        # assuming the sub-* naming convention is used
        nm = ('_' + str(str(list(filter(lambda x: ('sub') in x, (in1.split('/'))))[0]).split('_')[0]).split('-')[1])

        # handle workflows
        if wf == 1:
            # to get the count, indices and coordinates of overlapping voxels
            ov_count = np.count_nonzero(in_overlap)
            ov_idx = np.where(in_overlap)
            ov_ijk = np.vstack(ov_idx).T
            ov_xyz = nib.affines.apply_affine(aff2, ov_ijk) # this actually works

            ov_cog_ijk = ndimage.measurements.center_of_mass(in_overlap)
            ov_cog_xyz = nib.affines.apply_affine(aff2, ov_cog_ijk)

            # create empty arrays for distances
            ov_cog_2_maskAv_ds = np.zeros((ijk1.shape[0]), np.float32) # for ov COG to mask A voxels
            OVv_2_maskACOG_ds = np.zeros((ov_ijk.shape[0]), np.float32) # for ov voxels to COGA
            OVv_2_maskBCOG_ds = np.zeros((ov_ijk.shape[0]), np.float32) # for ov voxels to COGB

            # create empty arrays for voxel maps
            ov_cog_2_maskAv_map = np.zeros(im1_data.shape, np.uint16)
            OVv_2_maskACOG_map = np.zeros(in_overlap.shape, np.uint16)
            OVv_2_maskBCOG_map = np.zeros(in_overlap.shape, np.uint16)

            ov_Vox_map = np.zeros(in_overlap.shape, np.uint16) # for mapping the overlapping voxels to image
            ov_Vox_COG_map = np.zeros(in_overlap.shape, np.uint16) # for mapping COG of overlapping voxels


            # loop 1 to get dist. between all maskA voxels and ov_COG
            for ww in range(0,xyz1.shape[0]):
                ov_cog_2_maskAv_ds[ww] = np.linalg.norm(xyz1[ww]-ov_cog_xyz)


            # loop 2 to get dist. between all ov voxels and maskA_COG
            for ff in range(0,ov_xyz.shape[0]):
                OVv_2_maskACOG_ds[ff] = np.linalg.norm(ov_xyz[ff]-cog1_xyz)
                OVv_2_maskBCOG_ds[ff] = np.linalg.norm(ov_xyz[ff]-cog2_xyz)

            # find mins
            min_AvsOVCOG_d = np.amin(ov_cog_2_maskAv_ds)
            min_OV2ACOG_d = np.amin(OVv_2_maskACOG_ds)
            min_OV2BCOG_d = np.amin(OVv_2_maskBCOG_ds)

            # to get array indices of min values
            idx_1 = np.where(ov_cog_2_maskAv_ds == min_AvsOVCOG_d)
            idx_2 = np.where(OVv_2_maskACOG_ds == min_OV2ACOG_d)
            idx_3 = np.where(OVv_2_maskBCOG_ds == min_OV2BCOG_d)

            # to get actual coordinates
            ca_vox_mm = xyz1[idx_1[0]]
            ca_vox_vv = ijk1[idx_1[0]]
            cb_vox_mm = ov_xyz[idx_2[0]]
            cb_vox_vv = ov_ijk[idx_2[0]]
            cc_vox_mm = ov_xyz[idx_3[0]]
            cc_vox_vv = ov_ijk[idx_3[0]]

            # create array for ov_COG 2 maskA_voxels
            ov_cog_2_maskAv_map[np.int16(ca_vox_vv[0][0]), np.int16(ca_vox_vv[0][1]), np.int16(ca_vox_vv[0][2])] = 1
            dil_11 = np.uint16(ndimage.morphology.binary_dilation(ov_cog_2_maskAv_map, iterations=5))
            dil_11[np.int16(ca_vox_vv[0][0]), np.int16(ca_vox_vv[0][1]), np.int16(ca_vox_vv[0][2])] = 10
            # create array for ov voxels to maskA COG
            OVv_2_maskACOG_map[np.int16(cb_vox_vv[0][0]), np.int16(cb_vox_vv[0][1]), np.int16(cb_vox_vv[0][2])] = 1
            dil_22 = np.uint16(ndimage.morphology.binary_dilation(OVv_2_maskACOG_map, iterations=5))
            dil_22[np.int16(cb_vox_vv[0][0]), np.int16(cb_vox_vv[0][1]), np.int16(cb_vox_vv[0][2])] = 10
            # create array for ov voxels to maskB COG
            OVv_2_maskBCOG_map[np.int16(cc_vox_vv[0][0]), np.int16(cc_vox_vv[0][1]), np.int16(cc_vox_vv[0][2])] = 1
            dil_33 = np.uint16(ndimage.morphology.binary_dilation(OVv_2_maskBCOG_map, iterations=5))
            dil_33[np.int16(cc_vox_vv[0][0]), np.int16(cc_vox_vv[0][1]), np.int16(cc_vox_vv[0][2])] = 10

            # Create array for ov_COG image and dilate
            ov_Vox_COG_map[np.int16(ov_cog_ijk[0]), np.int16(ov_cog_ijk[1]), np.int16(ov_cog_ijk[2])] = 1
            dilated_cogOV = np.uint16(ndimage.morphology.binary_dilation(ov_Vox_COG_map, iterations=5))
            dilated_cogOV[np.int16(ov_cog_ijk[0]), np.int16(ov_cog_ijk[1]), np.int16(ov_cog_ijk[2])] = 10

            # Create array for all ov_voxels image (will need a for loop)
            for qq in range(0,ov_xyz.shape[0]):
                ov_Vox_map[np.int16(ov_ijk[qq][0]), np.int16(ov_ijk[qq][1]), np.int16(ov_ijk[qq][2])] = 1

            nib.save(nib.Nifti1Image(np.uint16(ov_Vox_map), aff2), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_initial_overlapping_voxels.nii.gz')
            nib.save(nib.Nifti1Image(np.uint16(dilated_cogOV), aff2), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_initial_overlapping_voxels_COG.nii.gz')
            nib.save(nib.Nifti1Image(np.uint16(dil_11), aff1), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_maskA_vox_mindist_2_overlap_COG.nii.gz')
            nib.save(nib.Nifti1Image(np.uint16(dil_22), aff2), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_overlap_vox_mindist_2_mask_A_COG.nii.gz')
            nib.save(nib.Nifti1Image(np.uint16(dil_33), aff2), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_overlap_vox_mindist_2_mask_B_COG.nii.gz')

            # calc percent overlap
            ov_perc_mB = 100.0 * np.float32(np.count_nonzero(ov_ijk)) / np.float32(np.count_nonzero(im2_data))
            ov_perc_mA = 100.0 * np.float32(np.count_nonzero(ov_ijk)) / np.float32(np.count_nonzero(im1_data))

            print('Minimum distance between COG of overlapping voxels and all voxels of mask A = ', min_AvsOVCOG_d , 'mm')
            print('Minimum distance between COG of mask A and all overlapping voxels = ', min_OV2ACOG_d , 'mm')
            print('Minimum distance between COG of mask B and all overlapping voxels = ', min_OV2BCOG_d , 'mm')
            print('Number of overlapping voxels between both masks: ', str(ov_count), ' voxels')
            print('Percent volume overlap between masks relative to mask A = ', str(ov_perc_mA).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", ""), '%')
            print('Percent volume overlap between masks relative to mask B = ', str(ov_perc_mB).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", ""), '%')


            # need to propagate to results text file then append later results to it without overwriting
            # also need to save overlap to nifti image ;)
            with open(pwd + '/' + out_n + '_output/' +  out_n + '_output' + nm + '_output_measures.txt', "w" ) as file_handler:
                file_handler.write('Initial overlap found between both masks, distance calculations using overlapping voxels, their COG, as well as external outlines and COGs of both masks' + '\n' + \
                '\n' + \
                'Minimum distance between COG of overlapping voxels and all voxels of mask A = ' + str(min_AvsOVCOG_d).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + 'mm' + '\n' + \
                'Minimum distance between COG of mask A and all overlapping voxels = ' + str(min_OV2ACOG_d).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + 'mm' + '\n' + \
                'Minimum distance between COG of mask B and all overlapping voxels = ' + str(min_OV2BCOG_d).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + 'mm' + '\n' + \
                '\n' + \
                'Number of overlapping voxels between both masks: ' + str(ov_count) + ' voxels' + '\n' + \
                'Percent volume overlap between masks relative to mask A = ' + str(ov_perc_mA) + '%' + '\n' + \
                'Percent volume overlap between masks relative to mask B = ' + str(ov_perc_mB) + '%' + '\n' \
                '\n')
                file_handler.close()

        elif wf == 2:
            # print('No overlap found')
            # need to propagate to results text file then append later results to it without overwriting
            with open(pwd + '/' + out_n + '_output/' +  out_n + '_output' + nm + '_output_measures.txt', "w" ) as file_handler:
                file_handler.write('No overlap found between both masks, distance calculations done using external outlines and COGs only' + '\n' + '\n')
                file_handler.close()

        # Save intermediate images to nii.gz in output dir
        nib.save(nib.Nifti1Image(outline1, aff1), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_mask_A_edge.nii.gz')
        nib.save(nib.Nifti1Image(outline2, aff2), pwd +  '/' + out_n + '_output' + '/' + out_n + nm + '_mask_B_edge.nii.gz')

        # needs a better cleanup strategy than simple morpho closure
        # potential helpful option -> https://www.delftstack.com/howto/python/smooth-data-in-python/
        # nib.save(nib.Nifti1Image(clean_im1, aff1), pwd +  '/' + out_n + '_output' + '/' + out_n + nm + '_mask_A_cleaned.nii.gz')
        # nib.save(nib.Nifti1Image(clean_im2, aff2), pwd +  '/' + out_n + '_output' + '/' + out_n + nm + '_mask_B_cleaned.nii.gz')
        
        nib.save(nib.Nifti1Image(im1_data, aff1), pwd +  '/' + out_n + '_output' + '/' + out_n + nm + '_mask_A.nii.gz')
        nib.save(nib.Nifti1Image(im2_data, aff2), pwd +  '/' + out_n + '_output' + '/' + out_n + nm + '_mask_B.nii.gz')

        nib.save(nib.Nifti1Image(eroded_im1, aff1), pwd +  '/' + out_n + '_output' + '/' + out_n + nm + '_mask_A_eroded.nii.gz')
        nib.save(nib.Nifti1Image(eroded_im2, aff2), pwd +  '/' + out_n + '_output' + '/' + out_n + nm + '_mask_B_eroded.nii.gz')

        # save output measures to a text file
        with open(pwd + '/' + out_n + '_output/' +  out_n + '_output' + nm + '_output_measures.txt', "a+" ) as file_handler:
            file_handler.write('Minimum distance between all voxels of mask A and mask B: ' + \
                str(all_min) + 'mm \n' + \
                'This is found between:- ' + '\n' + \
                'Mask A voxel at voxel coordinates: ' + str(a_vox_vv).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                'Mask A voxel at mm coordinates: ' + str(a_vox_mm).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' \
                'Mask B voxel at voxel coordinates: ' + str(b_vox_vv).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                'Mask B voxel at mm coordinates: ' + str(b_vox_mm).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                '\n' + \
                'Minimum distance between COG of mask A and COG of mask B: ' + str(cogs_d) + 'mm \n' + \
                'COG of mask A voxel coordinates: ' + str(cog1).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                'COG of mask A mm coordinates: ' + str(cog1_xyz).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                'COG of mask B voxel coordinates :' + str(cog2).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                'COG of mask B mm coordinates: ' + str(cog2_xyz).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                '\n' + \
                'Minimum distance between COG of mask A and all voxels of mask B: ' + str(coga_2b) + 'mm \n' + \
                'Mask B voxel(s) with shortest distance to mask A COG voxel coordinates: ' + str(ca2bijk).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                'Mask B voxel(s) with shortest distance to mask A COG mm coordinates: ' + str(ca2bxyz).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                '\n' + \
                'Minimum distance between COG of mask B and all voxels of mask A: ' + str(cogb_2a) + 'mm \n' + \
                'Mask A voxel(s) with shortest distance to mask B COG voxel coordinates: ' + str(cb2aijk).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n' + \
                'Mask A voxel(s) with shortest distance to mask B COG mm coordinates: ' + str(cb2axyz).replace("  ", " ").replace(" ", ", ").replace("[", "").replace("]", "") + '\n')
                
                
        # save voxels of min distances to two different images
        # save voxels of min distances to the same image or different images ??
        vox_A_maps[np.int16(a_vox_vv[0][0]), np.int16(a_vox_vv[0][1]), np.int16(a_vox_vv[0][2])] = 1
        COGA_map[np.int16(cog1[0]), np.int16(cog1[1]), np.int16(cog1[2])] = 1
        vox_B_maps[np.int16(b_vox_vv[0][0]), np.int16(b_vox_vv[0][1]), np.int16(b_vox_vv[0][2])] = 1
        COGB_map[np.int16(cog2[0]), np.int16(cog2[1]), np.int16(cog2[2])] = 1
        dilated_A = np.uint16(ndimage.morphology.binary_dilation(vox_A_maps, iterations=5))
        dilated_cogA = np.uint16(ndimage.morphology.binary_dilation(COGA_map, iterations=5))
        dilated_B = np.uint16(ndimage.morphology.binary_dilation(vox_B_maps, iterations=5))
        dilated_cogB = np.uint16(ndimage.morphology.binary_dilation(COGB_map, iterations=5))

        dilated_A[np.int16(a_vox_vv[0][0]), np.int16(a_vox_vv[0][1]), np.int16(a_vox_vv[0][2])] = 10
        dilated_B[np.int16(b_vox_vv[0][0]), np.int16(b_vox_vv[0][1]), np.int16(b_vox_vv[0][2])] = 10
        dilated_cogA[np.int16(cog1[0]), np.int16(cog1[1]), np.int16(cog1[2])] = 10
        dilated_cogB[np.int16(cog2[0]), np.int16(cog2[1]), np.int16(cog2[2])] = 10
        
        # save these voxel maps
        nib.save(nib.Nifti1Image(np.uint16(dilated_A), aff1), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_mask_A_vox_mindist_2_all_B_mask_vox.nii.gz')
        nib.save(nib.Nifti1Image(np.uint16(dilated_B), aff2), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_mask_B_vox_mindist_2_all_A_mask_vox.nii.gz')
        nib.save(nib.Nifti1Image(np.uint16(dilated_cogA), aff1), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_mask_A_COG.nii.gz')
        nib.save(nib.Nifti1Image(np.uint16(dilated_cogB), aff2), pwd + '/' + out_n + '_output' + '/' + out_n + nm + '_mask_B_COG.nii.gz')

    else:         
        print('the affines of the inputs are not matching, please double check, exiting')
        exit()


if __name__ == "__main__":
   main(sys.argv[1:])
