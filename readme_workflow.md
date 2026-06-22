## Step 1 : convert dicom 

```bash
# Convert CT DICOM directory
radot dicom ct_dcm -o ct.nii.gz

# Convert SPECT DICOM file
radot dicom spect.dcm -o spect.nii.gz
```

## Step 2 : add/check metadata

```bash
radot meta ct.nii.gz --modality ct

radot meta spect.nii.gz --modality spect --weight 80 --activity 120 --injection-time "2026-05-22 10:30:00" --acquisition-time "2026-05-22 10:45:00"
```

## Step 3 : segment (and stats)

```bash
# Segment CT using TotalSegmentator with dynamic AAPM TG-263 remapping
radot segment ct.nii.gz -o labels.nii.gz

# Check with vvv
vvv ct.nii.gz + labels.nii.gz

# Compute stats on CT
radot stats labels.nii.gz -i ct.nii.gz

# Compute stats on SPECT
radot stats labels.nii.gz -i spect.nii.gz
```

## Step 4 : dose

Coming soon!