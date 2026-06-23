## Step 1 : convert dicom 

```bash

# list dicoms in some folders
radot dicom-list folder/

# Convert CT DICOM directory
radot dicom ct_dcm -o ct.nii.gz

# Convert SPECT DICOM file
radot dicom spect.dcm -o spect.nii.gz
```

## Step 2 : add/check metadata

```bash
# Print metadata
radot meta ct.nii.gz

# Set metadata
radot meta spect.nii.gz --modality spect --weight 98 --activity 7400 --injection-time "2021-10-20 09:00:00" --acquisition-time "2021-10-21 12:10:09"

radot stats -i spect.nii.gz --unit bq

```

## Step 3 : segment (and stats)

```bash
# Segment CT using TotalSegmentator with dynamic AAPM TG-263 remapping
radot segment ct.nii.gz -o labels.nii.gz

# Check with vvv
vvv ct.nii.gz + labels.nii.gz
vvv ct.nii.gz, spect.nii.gz, jet

# Compute stats on CT
radot stats labels.nii.gz -i ct.nii.gz

# Compute stats on SPECT
radot stats labels.nii.gz -i spect.nii.gz --unit bq

```

## Step 4 : dose

Coming soon!


(TODO later -> same workflow without meta, explicit convert)