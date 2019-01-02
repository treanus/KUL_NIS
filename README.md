# KUL_neuroimaging_tools
Tools to analyse (resting/task-based) fMRI, diffusion & structural MRI data

Project PI's: Stefan Sunaert

Contributors: A. Radwan

# testing_contributors

Requires Mrtrix3, FSL, ants, dcm2niix, dcm2bids (jooh fork), docker, fmriprep and mriqc

@ Stefan Sunaert - UZ/KUL - stefan.sunaert@uzleuven.be

v0.2a - dd 12/12/2018 - pre beta version

These tools are only intended for MRI processing at Stefan's lab

They will procces an entire study (multiple subjects) with structural, functional and diffusion data:

        - convert dicom files to BIDS format
        - perform mriqc on structural and functional data
        - perform fmriprep on structural and functional data
        - perform freesurfer on the structural data (only T1w for now) 
        - perform mrtix3 processing on dMRI data
        - optionally:
            - perform combined structural and dMRI data analysis (depends on fmriprep, freesurfer and mrtrix3 above)
            - perfrom dbsdrt (automated tractography of the dentato-rubro-thalamic tract) on dMRI + structural data (depends on all above)
 
 Requirements:
 
        A correct installation of your mac (for now, maybe later also a hpc) at the lab
            - including:
                - dcm2niix (in KUL_apps)
                - dcm2bids (jooh fork, using pip)
                - docker
                - freesurfer (in KUL_apps)
                - mrtrix (in KUL_apps)  
                - last but not least, a correct installation of up-to-date KUL_NeuroImaging_Tools (in KUL_apps)
                - correct setup of your .bashrc and .bash_profile
 
 It uses a major config file, e.g. "study_config/subjects_and_options.csv" in which one informs the tools about:
 
        What and how (options) to perform:
                - mriqc (yes/no) 
                    (no options implemented yet)
                - fmriprep (yes/no), and specifies options:
                    all fmriprep options may be given,
                    e.g.:
                        --anat-only (to only process structural)
                        --longitudinal (to process longitudinal data)
                        --ignore fieldmaps
                - freesurfer (yes/no) 
                    (no options implemented yet)
                - KUL_dwiprep processing, i.e. a full mrtrix processing pipeline (yes/no)
                    options may be e.g.:
                        --slm=linear --repol (to provide to eddy)
                - KUL_dwiprep_anat processing (yes/no) 
                    (no options implemented yet)
                - KUL_dwiprep_dbsdrt processing (yes/no)
                        option nods e.g. 4000
