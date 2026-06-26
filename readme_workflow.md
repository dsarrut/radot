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

## Step 3 bis (optional) ; resample if needed

```bash
radot resample ct.nii.gz -s 4 4 4 -o ct_resampled.nii.gz --unit HU

radot resample spect.nii.gz -l labels_all.nii.gz -o spect_resampled.nii.gz -t spect --unit suv
```

## Step 3 ter (optional) ; stats

```bash
radot stats labels_all.nii.gz -i spect.nii.gz --unit bqml --rois Liver,Kidney_L
```

## Step 4 : dose

Compute the dose with the command line: 

```bash
radot dose spect.nii.gz -l labels_all.nii.gz -m hanscheid2017  --isotope lu177 --resample like_activity --verbose --ct ct.nii.gz --teff Liver:45 --teff kidney_l:32 --rois Liver --rois Kidney_L
```

Or with the equivalent Python function: 

```python
results = compute_roi_dose_file(
    "spect.nii.gz", 
    "labels_all.nii.gz", 
    method="hanscheid2017",
    rois=["Liver", "Kidney_L"],
    isotope="Lu177",
    ct_path="ct.nii.gz",
    resample_target="like_activity",
    method_params={"Liver": {"teff": 45.0}, "Kidney_L": {"teff": 32.0}},
    verbose=True,
)
```

