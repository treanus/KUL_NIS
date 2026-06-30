# KUL_dcm2bids

## Purpose

Use this command to convert dicom files to the BIDS format.



## Organise your data first

Although KUL_dcm2bids works however you organise the data, we suggest to use the following folder/data organisation (example):

![Image](KUL_dcm2bids_1.png)

Note that the MCI_009 folder could have been compressed into 1 zip or tar.gz file.   
Note too that the dicom folder could have contained thousands of dicom files in an unorganised manner.


## Usage

A typical command for conversion could be:  

`KUL_dcm2bids.sh -p MCI009 -s 1 -d DICOM/MCI_009/ -e -v -c study_config/sequences.txt`



You need to be in **main** directory to run this command (in the example /DATA/My_Study). This is also where the BIDS folder will be created. Note that we did ***not*** use an underscore in the participant name (MCI009 not MCI_009, which is not BIDS compliant).

For info about the command just run:

`KUL_dcm2bids.sh`


And this will output:

```
...more info above

Example:

  KUL_dcm2bids.sh -p pat001 -d pat001.zip -c definitions_of_sequences.txt

Required arguments:

     -d:  dicom_zip_file (the zip or tar.gz containing all your dicoms, or directory containing dicoms)
     -p:  participant (anonymised name of the subject in bids convention)
     -c:  definitions of sequences (T1w=MPRAGE,dwi=seq, etc..., see above)

Optional arguments:

     -o:  bids directory
     -s:  session (for longitudinal study with multiple timepoints)
     -t:  temporary directory (default = /tmp)
     -e:  copy task-*_events.tsv from config to BIDS dir
     -a:  further anonymise the subject by using pydeface (takes much longer)
     -v:  verbose 
```

Note: you can  specify a directory containing dicoms, but also use a zip or tar.gz file. The latter can be more handy (especially if you have to transfer thousands of them from/to a USB stick e.g.).

When you converted data (one or multiple subjects and sessions) please use the bids validator to check your BIDS directory.

## Output

Following the example above you will get:

![Image](KUL_dcm2bids_2.png)


## Script variants (dcm2bids v3 vs v2)

The dcm2bids tool changed its configuration schema between major versions. KUL_NIS ships three variants so you can match whichever dcm2bids is installed on your system:

| Script | dcm2bids version | Config schema | When to use |
|---|---|---|---|
| **`KUL_dcm2bids.sh`** (default) | **v3 (>= 3.0)** | `dataType` / `modalityLabel` / `customLabels` / `sidecarChanges` | Standard, recommended. This is the maintained version. |
| `KUL_dcm2bids_dev.sh` | v3 (>= 3.0) | same as above, **plus Philips PE-direction auto-detection** from DICOM private tag `[2005,107B]` | Testing only. Auto-derives phase-encoding direction for Philips DWI when header info is present. Not yet validated for all cases. |
| `KUL_dcm2bids_v2bkup.sh` | v2 (< 3.0) | `custom_entities` / `sidecar_changes` | Fallback for sites still pinned to dcm2bids v2. |

`KUL_dcm2bids.sh` runs a version-detection check at startup: if it finds dcm2bids v2 installed it prints a warning and points you to `KUL_dcm2bids_v2bkup.sh`. Upgrade with `pip install --upgrade dcm2bids`.

**Your config file schema must match the installed dcm2bids version** — a v2 config will not work with v3 and vice versa.

## Config file

The KUL_dcm2bids config file describes what data you have acquired on the scanner and how it should be converted into the BIDS format.

A typical config file could look as follows:

```
#Identifier,search-string,task,mb,pe_dir,acq_label  
T1w,MPRAGE  
cT1w,T1_Gd  
FLAIR,3D_FLAIR_mind_study  
func,rsfMRI,rest,8,j,multiTE  
func,tb_fMRI,nback,2,j  
dwi,part1,-,2,j,b2000  
dwi,part2,-,2,j,b4000  
dwi,revphase,-,2,j-,rev
```  


The first column describes the BIDS data type. This is the type of scan that you have acquired, e.g. a "FLAIR", or a "func" for functional MRI data.
In the example above, we acquired a T1w (without intravenous Gd contrast), a cT1w (with contrast) and a FLAIR image, two functional MRI scans (1 multi echo resting-state and 1 singel echo during an n-back task) and the DTI data were acquired using 3 diffusion sequences. In a first diffusion acquisition b0 and b2000 were acquired. In the second diffusion scan b0 and b4000 were acquired. In a third scan a few b0 and 6 directions with b2000 weigthing were again acquired, but now using "reverse phase-encoding). 



KUL_dcm2bids needs information how to identify e.g. the diffusion scan with b0 and b4000 data. Therefore we use the second column. This second column defines a search pattern, i.e. a unique string that will search in the Series Description for e.g. the "FLAIR" or functional MRI data that you acquired. The Series Description was defined on the scanner (and corresponds to a dicom tag). On Philips scanners it corresponds to what you defined in the ExamCard. It could e.g. be "Smartbrain" or "3D_MPRAGE" or "DTI_part1_optimised_dd_21032017". 
KUL_dcm2bids will use the second column (I.e. the search pattern) to search through all dicom files to find e.g. a FLAIR.
Note that the search-pattern needs to be unique to define a certain data type. 



The third column gives information about the fMRI task. 



The config file requires more information for the data type "func" acquired on Philips scanners:

- Column 4 = mb: defines the multiband factor used in the GE-EPI sequence
- Column 5 = pe_dir: defines the phase-encoding direction used in the EPI sequence:
  -  for phase-encoding direction AP, fat-shift P: use j
  -  for phase-encoding direction AP, fat-shift A: use j-

The config file requires more information for the data type "dwi" acquired on Philips scanners:

- Column 3 = task: obsolote here (set to -)
- Column 4 = mb: defines the multiband factor used in the SE-EPI sequence
- Column 5 = pe_dir: defines the phase-encoding direction used in the EPI sequence

Note that the config file ONLY needs the above Columns 4 and 5 on PHILIPS Scanners. Siemens and GE have this encoded in the dicom header.

An example config file can be found in [study_config/sequences.txt](/study_config/sequences.txt)


## Supported images

### Structural

- T1w, cT1w, FLAIR, T2w, PDw: fully BIDS compliant
- FGATIR:  Fast Gray Matter Acquisition T1 Inversion Recovery, this is not yet specified in the BIDS
- SWI: susceptibility weigthed images. Magnitude (SWI) and phase (SWIp) are split into separate series; not yet specified in the BIDS (listed in `.bidsignore`)
- DIR: double inversion recovery; not yet specified in the BIDS (listed in `.bidsignore`)
- MP2RAGE: stored with its INV1/INV2/UNI parts; not yet specified in the BIDS (listed in `.bidsignore`)
- MTI: magnetisation transfer contrast images (a pair of images), this is not yet specified in the BIDS

### fMRI

- func: both single and multi echo data
- events.tsv: put these in the study_config file, and use option -e
- sbref: single band reference 
- fmap: a field map

### DWI

- dwi: diffusion MRI

### ASL

- ASL: arterial spin labeling, images that have already been processed on the MR scanner (a rCVF map), this is not what is specified in the BIDS convention (those are native images)


## Dependencies

Internally KUL_dcm2bids uses:
- [dcm2bids](https://github.com/UNFmontreal/Dcm2Bids) — **v3 (>= 3.0)** for the default script; v2 only with `KUL_dcm2bids_v2bkup.sh`
- [dcm2niix](https://github.com/rordenlab/dcm2niix)
- [pydeface](https://github.com/poldracklab/pydeface) — only when using `-a` (further anonymisation)

These could be installed using [KUL_Linux_Installation](https://github.com/treanus/KUL_Linux_Installation)
