#!/usr/bin/env python

import argparse
import glob 
import os


parser = argparse.ArgumentParser(description="Run the longitudinal version of samseg",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-v", "--verbose", action="store_true", help="increase verbosity")
parser.add_argument("-F", "--flair", action="store_true", help="use the flair as well")
parser.add_argument("-n", "--ncpu", help="number of threads to use for samseg")
parser.add_argument("dest", help="Destination location")
args = parser.parse_args()
config = vars(args)
#print(config)

if args.ncpu is None:
    ncpu = 15
else:
    ncpu = args.ncpu

bidsdir = './BIDS'
outdir = args.dest
#print(outdir)


for root, dirs, files in os.walk(bidsdir):
    for dir in dirs:
        if 'sub-' in dir:
            searchdir = os.path.join(root, dir)
            #print(searchdir)
            #print(dir)
            os.makedirs(os.path.join(outdir,dir), exist_ok=True)

            for ImType in ["T1w"]:

                searchIms = searchdir + '/*/anat/*' + ImType + '.nii.gz'
                #print(searchIms)
                
                # find all Im
                Ims = glob.glob(searchIms)
                #print(Ims)

                Imsreg_input = [] 
                flairreg_input = []
                Imsreg_output = [] 
                flairreg_output = [] 
                for Im in Ims:
                    # regrid the input to 1mm isotropic
                    
                    dir_name, base_name = os.path.split(os.path.splitext(os.path.splitext(Im)[0])[0])
                    #print(dir_name)
                    base_name = base_name.split('_T1w')[0]
                    #print(base_name)
                    flair_search = os.path.join(dir_name, base_name + '_FLAIR.nii.gz')
                    T1w_iso_output = os.path.join(outdir, dir, base_name + '_T1w_iso.nii.gz')
                    #print(T1w_iso_output)
                    cmd = 'mrgrid ' + Im + ' regrid -voxel 1 ' + T1w_iso_output
                    print(cmd)
                    if not os.path.exists(T1w_iso_output):
                        out = os.popen(cmd).read().strip()
                        print(out)
                    #print(flair_search)
                    if args.flair: 
                        if os.path.exists(flair_search):
                            #print('There is a flair')
                            flair_iso_output = os.path.join(outdir, dir, base_name + '_FLAIR_iso.nii.gz')
                            #print(flair_iso_output)
                            cmd = 'mrgrid ' + flair_search + ' regrid -voxel 1 ' + flair_iso_output
                            #print(cmd)
                            if not os.path.exists(flair_iso_output):
                                out = os.popen(cmd).read().strip()
                                print(out)
                            Imsreg_input.append(T1w_iso_output)
                            flairreg_input.append(flair_iso_output)
                            Imsreg_output.append(os.path.join(outdir, dir, base_name + '_T1w_iso_reg.mgz'))
                            flairreg_output.append(os.path.join(outdir, dir, base_name + '_FLAIR_iso_reg.mgz'))
                    else:
                        Imsreg_input.append(T1w_iso_output)
                        Imsreg_output.append(os.path.join(outdir, dir, base_name + '_T1w_iso_reg.mgz'))
                #print(Imsreg_input)
                #print(Imsreg_output)
                
                mean_template = os.path.join(outdir, dir, 'T1w_mean.mgz')
                if not os.path.exists(mean_template):
                    cmd = 'mri_robust_template --mov ' + ' '.join(Imsreg_input) + ' --template ' + mean_template \
                        + ' --satit --mapmov ' + ' '.join(Imsreg_output)
                    print(cmd)

                    out = os.popen(cmd).read().strip()
                    print(out)

                if args.flair:
                    i=0
                    for Im in flairreg_input:
                        dir_name, base_name = os.path.split(os.path.splitext(os.path.splitext(Im)[0])[0])
                        lta = os.path.join(outdir, dir, base_name + '_FLAIRtoT1.lta')
                        if not os.path.exists(lta):
                            cmd = 'mri_coreg --threads ' + str(ncpu) + ' --mov ' + Im  + ' --ref ' + Imsreg_output[i] + ' --reg ' + lta
                            print(cmd)
                            out = os.popen(cmd).read().strip()
                            print(out)
                            cmd = 'mri_vol2vol --mov ' + Im + ' --reg ' + lta + ' --o ' + flairreg_output[i] + ' --targ ' + Imsreg_output[i]
                            print(cmd)
                            out = os.popen(cmd).read().strip()
                            print(out)
                        i = i + 1
                
                cmd_input = []
                i=0
                for samseg_input in Imsreg_output:
                    if args.flair:
                        cmd_input.append('--timepoint ' + flairreg_output[i] + ' ' + Imsreg_output[i])
                    else: 
                        cmd_input.append('--timepoint ' + samseg_input)
                    i = i + 1
                cmd = 'run_samseg_long ' + ' '.join(cmd_input) + ' --output ' + os.path.join(outdir, dir, 'samseg') + \
                    ' --lesion --lesion-mask-pattern 1 0 --threshold 0.7 ' + ' --threads ' + str(ncpu)
                if not os.path.exists(os.path.join(outdir, dir, 'samseg')):
                    print(cmd)
                    out = os.popen(cmd).read().strip()
                    print(out)
                
