#!/bin/bash
#
# Interactive helper for adding a single fMRI activation label to an existing
# Karawun folder by hand. NOT part of the automatic pipeline: KUL_clinical_fmridti.sh
# writes the fMRI labels itself (search for _fmri_label_colors) when run with -R.
# Use this only for one-off additions, and keep its palette values in sync with
# the reserved pool there.

extension=".nii"  # Replace with your desired extension
directory="."  # Replace with your desired directory

# Function to display the menu
display_menu() {
    local count=0
    for file in "${files[@]}"; do
        echo "$count. $file"
        ((count++))
    done
}

# Read files in the directory with the specified extension and store them in an array
files0=()
while IFS= read -r -d '' file; do
    files0+=("$file")
done < <(find "$directory" -type f -name "*$extension" -print0)

# Sort the files
IFS=$'\n' files=($(sort -n <<<"${files0[*]}"))

# Display the menu
echo "Files with extension '$extension' in the directory:"
display_menu

# Prompt for user input (file selection)
read -rp "Enter the number of the file you want to select: " choice

# Validate the file choice
if [[ $choice =~ ^[0-9]+$ ]] && ((choice >= 0 && choice < ${#files[@]})); then
    selected_file="${files[$choice]}"
    echo "You selected: $selected_file"
else
    echo "Invalid file choice."
    exit 1
fi

# Prompt for user input (category selection)
categories=("afMRI_HAND" "afMRI_LIP" "afMRI_FOOT" "afMRI_TAAL" "rsfMRI_HAND" "rsfMRI_LIP" "rsfMRI_FOOT" "rsfMRI_TAAL")
echo "Categories:"
for ((i=0; i<${#categories[@]}; i++)); do
    echo "$i. ${categories[$i]}"
done

read -rp "Enter the number of the category: " category_choice

# Validate the category choice
if [[ $category_choice =~ ^[0-9]+$ ]] && ((category_choice >= 0 && category_choice < ${#categories[@]})); then
    selected_category="${categories[$category_choice]}"
    echo "Selected category: $selected_category"
else
    echo "Invalid category choice."
    exit 1
fi

# Set values based on selected category.
#
# These are Brainlab palette indices, and they must stay out of the ranges the
# rest of the pipeline uses (see KUL_karawun_prepare.sh for the full budget):
#   1-41 known tracts, 42-49 auto tracts, 50 lesion, 51-63 fMRI.
#
# They used to be 20-27, every one of which collided: 20/21/22/25/26/27 are
# tract colours and 23/24 are the DBS STN VOIs. So a hand-made afMRI_TAAL label
# came out the same colour as the left STN VOI, and afMRI_HAND the same as a
# tract. Now drawn from the same reserved pool KUL_clinical_fmridti.sh uses,
# ordered by measured CIEDE2000 separation from the tract colours.
#
# All of these are above 30, so they REQUIRE the extended-palette karawun fork;
# stock karawun clamps anything above 30 to a single entry.
case $selected_category in
    "afMRI_HAND")
        value=51
        ;;
    "afMRI_LIP")
        value=52
        ;;
    "afMRI_FOOT")
        value=53
        ;;
    "afMRI_TAAL")
        value=54
        ;;
    "rsfMRI_HAND")
        value=55
        ;;
    "rsfMRI_LIP")
        value=56
        ;;
    "rsfMRI_FOOT")
        value=57
        ;;
    "rsfMRI_TAAL")
        value=58
        ;;
    *)
        echo "Invalid category."
        exit 1
        ;;
esac

# Prompt for user input (threshold)
read -rp "Enter the threshold: " threshold

# Perform further actions based on the selected file, category, value, and threshold
# Add your code here to process the selected file, category, value, and threshold as needed
echo "File: $selected_file"
echo "Category: $selected_category"
echo "Value: $value"
echo "Threshold: $threshold"


cmd1="mrcalc $selected_file $threshold -gt $value -mul \$0/${selected_category}_thres_${threshold}.nii.gz"
cmd2="find ../../../Karawun/ -type d -name "labels" -exec bash -c '$cmd1' {} \;"
eval $cmd2

