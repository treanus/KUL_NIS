#!/usr/bin/env python
# Convert dicom to nifti using mrconvert
# Stefan Sunaert - 17/05/2023

import os
import argparse

# Get commandline
parser = argparse.ArgumentParser(description="Convert dicom to nifti",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--participant', help='participant id', required=True)
parser.add_argument('--dicomdir', help='dicom input directory', required=True)
parser.add_argument('--seriesnumbers', nargs='+', help='series numbers', required=True)
parser.add_argument('--type', nargs='+', help='type, e.g T1w, dwi, func', required=True)
parser.add_argument('--donor_dcm', nargs='+', help='the donor dicom used to extract tags, e.g. IM-001.dcm', required=True)
parser.add_argument('-i', '--inputtype', nargs='+', help='inputtype, e.g. M_FFE')
parser.add_argument('-o', '--outputtype', nargs='+', help='outputtype, e.g. phase')
parser.add_argument('-a', '--acquisition', nargs='+', help='acquisition, e.g. ap')
parser.add_argument('-e', '--pe_direction', nargs='+', help='phase encoding direction, e.g. j-')

args = parser.parse_args()

# set inputs and check
participant = args.participant
dcm_dir = args.dicomdir
dcm_series = args.seriesnumbers
bids_type = args.type
dcm_types = args.inputtype
nii_parts = args.outputtype
dcm_pe = args.acquisition
nii_pe = args.pe_direction
donor_dcm = args.donor_dcm

'''
print(participant)
print(dcm_dir)
print(dcm_series)
print(dcm_types)
print(bids_type)
print(nii_parts)
'''

# tags to get from donor
# a function to get tags
def getDicomTag(tag):
    cmd = 'dcminfo -tag ' + tag + ' ' + '\"' + donor_dcm[0] +'\"'
    #print(cmd)
    out = os.popen(cmd).read().strip()
    #print(out)
    return out.split(' ')[1]

# define tags to read from the donor dcm
dict_tags = {'Modality': '0008 0060', \
        'MagneticFieldStrength': '0018 0087', \
        'ImagingFrequency': '0018 0084', \
        'Manufacturer': '0008 0070', \
        'InstitutionName': '0008 0080', \
        'InstitutionAddress': '0008 0081', \
        'InstitutionalDepartmentName': '0008 1040', \
        'DeviceSerialNumber': '0018 1000', \
        'StationName': '0008 1010', \
        'BodyPartExamined': '0018 0015', \
        'PatientPosition': '0018 5100', \
        'SoftwareVersions': '0018 1020', \
        'MRAcquisitionType': '0018 0023', \
        'SeriesDescription': '0008 103E', \
        'ProtocolName': '0018 1030', \
        'WaterFatShift': '2001 1022', \
        'EPIFactor': '2001 1013', \
        'Rows': '0028 0010'}

# get the relevant tags
dict_dcm = {}
for key in dict_tags:
    print(key)
    print(dict_tags[key])
    dict_tags[key] = getDicomTag(dict_tags[key])
    #print(dict_tags[key])
    dict_dcm.update({key : dict_tags[key]})
#print(dict_dcm)

# calculate 
#ActualEchoSpacing = WaterFatShift / (ImagingFrequency * 3.4 * (EPI_Factor + 1))
#TotalReadoutTIme = ActualEchoSpacing * EPI_Factor
# EffectiveEchoSpacing = TotalReadoutTime / (ReconMatrixPE - 1)
ActualEchoSpacing = float(dict_dcm['WaterFatShift']) \
    / (float(dict_dcm['ImagingFrequency']) * 3.4 * (float(dict_dcm['EPIFactor']) + 1))
TotalReadoutTime = ActualEchoSpacing * float(dict_dcm['EPIFactor'])
EffectiveEchoSpacing = TotalReadoutTime / (float(dict_dcm['Rows']) - 1 )

'''
print(ActualEchoSpacing)
print(TotalReadoutTime)
print(EffectiveEchoSpacing)
'''

# insert into dict
dict_dcm.update({'TotalReadoutTime': TotalReadoutTime})
dict_dcm.update({'EffectiveEchoSpacing': EffectiveEchoSpacing})
print(dict_dcm)

# make the addition properties to insert to the mif or nii
additional_properties = ''
for key in dict_dcm:
    print(key)
    additional_properties = additional_properties + '-set_property ' + key + ' ' + str(dict_dcm[key]) + ' '
print(additional_properties)


#dcm_types = ['M_SE','I_SE','R_SE','PHASE']
#nii_parts = ['mag','imag','real','phase']
#dcm_types = ['M_SE']
#nii_parts = ['mag']
#dcm_pe = ['ap', 'pa']
#nii_pe = ['j-', 'j']

if not os.path.exists(dcm_dir):
    print(dcm_dir + ' does not exist')
    exit(1)

#nii_dir = os.path.join('BIDS/sub-' + participant, 'dwi')
nii_dir = '.'
if not os.path.exists(nii_dir):
   os.makedirs(nii_dir)

#for i, dcm_serie in enumerate(dcm_series):
i = 0
dcm_serie = dcm_series
#print(i)
for j, nii_part in enumerate(nii_parts):
    #j = 0
    #nii_part = nii_parts[0]
    #print(j)
    print(bids_type[i])
    if bids_type[i] == 'dwi':
        part = '_part-' + str(nii_part)
    else:
        part = ''

    if args.pe_direction:
        pe = '_acq-' + dcm_pe[i]
    else:
        pe = ''
    nii_file = 'sub-' + participant + pe + part + '_' + str(bids_type[i])
    nii_nii = os.path.join(nii_dir,nii_file) + '.nii.gz'
    nii_json = os.path.join(nii_dir,nii_file) + '.json'
    if bids_type[i] == 'dwi':
        nii_bval = os.path.join(nii_dir,nii_file) + '.bval'
        nii_bvec = os.path.join(nii_dir,nii_file) + '.bvec'
        export_grad_fls = " -export_grad_fsl " + nii_bvec + " " + nii_bval 
        property_pe = " -set_property \"PhaseEncodingDirection\" \"" + str(nii_pe[i]) + "\" "
    else:
        export_grad_fls = ""
        property_pe = ""

    # check if we need to check on a type (e.g. M_FFE) too & set the output type if not specified
    if dcm_types:
        search_type = " | grep " + dcm_types[j]    
    else:
        search_type = ""

    cmd = "echo q | mrinfo \"" + str(dcm_dir) + "\" 2>&1 | grep " + str(dcm_serie[0]) + \
        search_type + " | awk '{print $1}' | mrconvert " + \
        " \"" + str(dcm_dir) + "\" " +  \
        additional_properties + ' ' + \
        " -json_export " + nii_json + ' ' + \
        export_grad_fls + ' ' + \
        property_pe + ' ' + \
        nii_nii + ' -force'
    print(cmd)
    #exit()
    out = os.popen(cmd).read().strip()
    print(out)
