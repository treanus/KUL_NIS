#!/bin/bash

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
files=()
while IFS= read -r -d '' file; do
    files+=("$file")
done < <(find "$directory" -type f -name "*$extension" -print0)

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

# Set values based on selected category
case $selected_category in
    "afMRI_HAND")
        value=20
        ;;
    "afMRI_LIP")
        value=21
        ;;
    "afMRI_FOOT")
        value=22
        ;;
    "afMRI_TAAL")
        value=23
        ;;
    "rsfMRI_HAND")
        value=24
        ;;
    "rsfMRI_LIP")
        value=25
        ;;
    "rsfMRI_FOOT")
        value=26
        ;;
    "rsfMRI_TAAL")
        value=27
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

