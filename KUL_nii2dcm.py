#!/usr/bin/env python

# Convert nifti to dicom given a donor dicom image
# Stefan Sunaert - 27/02/2023
# Mainly based on SimpleITK - https://simpleitk.readthedocs.io/en/master/link_DicomSeriesFromArray_docs.html

import SimpleITK as sitk
import argparse
import sys
import time
import os
import shutil
import numpy as np


# Get and check commandline
parser = argparse.ArgumentParser(description="Convert a nifti or 3d-tiff to dicom given a donor dicom image",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-v", "--verbose", action="store_true", help="increase verbosity")
parser.add_argument("-s", "--seriesdescription")
parser.add_argument("nifti", help="nifti or 3d-tiff image")
parser.add_argument("donor", help="dicom donor image")
parser.add_argument("dicomdir", help="dicom output directory")
args = parser.parse_args()
config = vars(args)
#print(config)


# Define functions
def writeSlices(series_tag_values, new_img, out_dir, i):
    image_slice = new_img[:, :, i]

    # Tags shared by the series.
    list(
        map(
            lambda tag_value: image_slice.SetMetaData(
                tag_value[0], tag_value[1]
            ),
            series_tag_values,
        )
    )

    # Slice specific tags.
    #   Instance Creation Date
    image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))
    #   Instance Creation Time
    image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))

    # Setting the type to CT so that the slice location is preserved and
    # the thickness is carried over.
    #image_slice.SetMetaData("0008|0060", "MR")

    # (0020, 0032) image position patient determines the 3D spacing between
    # slices.
    #   Image Position (Patient)
    image_slice.SetMetaData(
        "0020|0032",
        "\\".join(map(str, new_img.TransformIndexToPhysicalPoint((0, 0, i)))),
    )
    #   Instance Number
    image_slice.SetMetaData("0020|0013", str(i))

    # Write to the output directory and add the extension dcm, to force
    # writing in DICOM format.
    writer.SetFileName(os.path.join(out_dir, str(i) + ".dcm"))
    writer.Execute(image_slice)



# set inputs and check
donor_dcm = args.donor
if not os.path.exists(donor_dcm):
    print(donor_dcm + ' does not exist')
    exit(1)
nifti_input = args.nifti
if not os.path.exists(nifti_input):
    print(nifti_input + ' does not exist')
    exit(1)
dcm_output = args.dicomdir
if args.seriesdescription:
    seriesdesc = args.seriesdescription
else:
    seriesdesc = 'IKTsimple - KUL_NIS'

# Read the donor DICOM
reader = sitk.ImageFileReader()
reader.SetFileName(donor_dcm)
reader.LoadPrivateTagsOn()
reader.ReadImageInformation()

# Display the tags
if args.verbose:
    for k in reader.GetMetaDataKeys():
        v = reader.GetMetaData(k)
        try: 
            print(f'({k}) = = "{v}"')
        except:
            print("An exception occurred")

# Copy relevant tags from the original meta-data dictionary (private tags are
# also accessible).
tags_to_copy = [
    "0010|0010",  # Patient Name
    "0010|0020",  # Patient ID
    "0010|0030",  # Patient Birth Date
    "0010|0040",  # Patient Sex
    "0020|000D",  # Study Instance UID, for machine consumption
    "0020|0010",  # Study ID, for human consumption
    "0008|0020",  # Study Date
    "0008|0030",  # Study Time
    "0008|0050",  # Accession Number
    "0008|0060",  # Modality
]

# Read the nii or tiff
nii_img = sitk.ReadImage(nifti_input)

# Convert the data to int16
np.img_data = sitk.GetArrayFromImage(nii_img)
max = np.amax(np.img_data)
#print(max)
img_int16 = np.img_data * ( np.iinfo(np.int16).max / max )
img_int16b = img_int16.astype(np.int16)
new_img = sitk.GetImageFromArray(img_int16b)
new_img.CopyInformation(nii_img)

'''
# Check the data type and set spacing in case of TIFF
try:
    print(nii_img.GetMetaData('nifti_type'))
    print('Input is a nifti')
except:
    print('Input is not nifti, probably TIFF; setting spacing to 1,1,1')
    new_img.SetSpacing([1.0, 1.0, 1.0])
'''

# Write the 3D image as a series
# IMPORTANT: There are many DICOM tags that need to be updated when you modify
#            an original image. This is a delicate operation and requires
#            knowledge of the DICOM standard. This example only modifies some.
#            For a more complete list of tags that need to be modified see:
#                  http://gdcm.sourceforge.net/wiki/index.php/Writing_DICOM
#            If it is critical for your work to generate valid DICOM files,
#            It is recommended to use David Clunie's Dicom3tools to validate
#            the files:
#                  http://www.dclunie.com/dicom3tools.html

writer = sitk.ImageFileWriter()
# Use the study/series/frame of reference information given in the meta-data
# dictionary and not the automatically generated information from the file IO
writer.KeepOriginalImageUIDOn()

modification_time = time.strftime("%H%M%S")
modification_date = time.strftime("%Y%m%d")

# Copy some of the tags and add the relevant tags indicating the change.
# For the series instance UID (0020|000e), each of the components is a number,
# cannot start with zero, and separated by a '.' We create a unique series ID
# using the date and time. Tags of interest:
direction = new_img.GetDirection()
series_tag_values = [
    (k, reader.GetMetaData(k))
    for k in tags_to_copy
    if reader.HasMetaDataKey(k)
] + [
    ("0008|0031", modification_time),  # Series Time
    ("0008|0021", modification_date),  # Series Date
    ("0008|0008", "DERIVED\\SECONDARY"),  # Image Type
    (
        "0020|000e",
        "1.2.826.0.1.3680043.2.1125."
        + modification_date
        + ".1"
        + modification_time,
    ),  # Series Instance UID
    (
        "0020|0037",
        "\\".join(
            map(
                str,
                (
                    direction[0],
                    direction[3],
                    direction[6],
                    direction[1],
                    direction[4],
                    direction[7],
                ),
            )
        ),
    ),  # Image Orientation
    # (Patient)
    ("0008|103e", seriesdesc),  # Series Description
]

# Give info
print('Incorporating the following dicom tags:')
print(series_tag_values)

# Ckean and Make the output dir
if os.path.exists(dcm_output):
    shutil.rmtree(dcm_output)
os.makedirs(dcm_output, exist_ok=True)

'''
if pixel_dtype == np.float64:
    # If we want to write floating point values, we need to use the rescale
    # slope, "0028|1053", to select the number of digits we want to keep. We
    # also need to specify additional pixel storage and representation
    # information.
    rescale_slope = 0.001  # keep three digits after the decimal point
    series_tag_values = series_tag_values + [
        ("0028|1053", str(rescale_slope)),  # rescale slope
        ("0028|1052", "0"),  # rescale intercept
        ("0028|0100", "16"),  # bits allocated
        ("0028|0101", "16"),  # bits stored
        ("0028|0102", "15"),  # high bit
        ("0028|0103", "1"),
    ]  # pixel representation
'''

# Write slices to output directory
list(
    map(
        lambda i: writeSlices(series_tag_values, new_img, dcm_output, i),
        range(new_img.GetDepth()),
    )
)

sys.exit(0)