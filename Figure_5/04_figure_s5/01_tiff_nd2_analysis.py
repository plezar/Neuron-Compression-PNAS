# process tiff and nd2 data for fosb analysis
# uses fixed thresholds from previous optimization

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from nd2reader import ND2Reader
from skimage import measure
from skimage.segmentation import clear_border
from skimage.morphology import remove_small_objects
import tifffile
from matplotlib.path import Path as mplPath
from scipy import stats
from datetime import datetime

# thresholds - these work well for our data
DAPI_THRESHOLD = 3363
GFAP_THRESHOLD = 900
MIN_SIZE = 300

def create_tiff_nd2_mapping():
    # map tiff slices to nd2 series
    nd2_files = {
        '176130_C': {'group': 'Control', 'fovs': 7},
        '178992_C': {'group': 'Control', 'fovs': 5},
        '178993_C': {'group': 'Control', 'fovs': 7},
        '179041_E': {'group': 'Experimental', 'fovs': 8},
        '179042_E': {'group': 'Experimental', 'fovs': 8}
    }
    
    mapping_data = []
    for replicate_id, info in nd2_files.items():
        for fov_num in range(1, info['fovs'] + 1):
            mapping_data.append({
                'replicate_id': replicate_id,
                'group': info['group'],
                'series_number': fov_num,
                'series_id': f"{replicate_id}_Series_{fov_num}",
                'analysis_strategy': 'TIFF_nuclei_template_applied_to_ND2_FosB',
                'nuclei_source': 'TIFF_MIP_with_exact_fixed_raw_thresholds',
                'fosb_source': 'ND2_all_z_slices'
            })
    
    return pd.DataFrame(mapping_data)

def get_biological_replicates_exact_method(input_dir):
    # find all tiff files and group by replicate
    replicates = {}
    
    for file in os.listdir(input_dir):
        if not file.endswith('.tif'):
            continue
            
        parts = file.split('_')
        replicate_id = parts[0]  # e.g. 176130
        
        if replicate_id not in replicates:
            replicates[replicate_id] = {'DAPI': None, 'FOSB': None, 'GFAP': None, 'NeuN': None}
        
        # figure out which channel this is
        file_upper = file.upper()
        if 'DAPI' in file_upper:
            replicates[replicate_id]['DAPI'] = os.path.join(input_dir, file)
        elif 'FOSB' in file_upper:
            replicates[replicate_id]['FOSB'] = os.path.join(input_dir, file)
        elif 'GFAP' in file_upper:
            replicates[replicate_id]['GFAP'] = os.path.join(input_dir, file)
        elif 'NEUN' in file_upper:
            replicates[replicate_id]['NeuN'] = os.path.join(input_dir, file)
    
    return replicates

def load_roi_exclusion_mask(replicate_id, slice_num, roi_dir):
    # load roi if it exists
    roi_json = os.path.join(roi_dir, f"{replicate_id}_slice_{slice_num}_exclusion_roi.json")
    
    if os.path.exists(roi_json):
        with open(roi_json, 'r') as f:
            pts = json.load(f)
        if pts and len(pts) > 2:
            return pts
    
    return None

def create_exclusion_mask(roi_points, image_shape):
    if roi_points is None:
        return np.zeros(image_shape, dtype=bool)
    
    y, x = np.mgrid[:image_shape[0], :image_shape[1]]
    points = np.vstack((x.ravel(), y.ravel())).T
    
    path = mplPath(roi_points)
    mask = path.contains_points(points)
    mask = mask.reshape(image_shape)
    
    return mask

def create_tiff_nuclei_template_exact_method(tiff_files, replicate_id, series_num, roi_dir, output_dir):
    # create template from tiff for this specific series
    print(f"  Creating TIFF nuclei template for {replicate_id} Series {series_num}...")
    
    # load tiff stacks
    dapi_stack = tifffile.imread(tiff_files['DAPI'])
    fosb_stack = tifffile.imread(tiff_files['FOSB'])
    gfap_stack = tifffile.imread(tiff_files['GFAP'])
    neun_stack = tifffile.imread(tiff_files['NeuN'])
    
    # get the right slice for this series
    slice_idx = series_num - 1
    
    if slice_idx >= dapi_stack.shape[0]:
        print(f"    ERROR: Series {series_num} not available in TIFF stack")
        return None
    
    dapi_slice = dapi_stack[slice_idx]
    fosb_slice = fosb_stack[slice_idx] 
    gfap_slice = gfap_stack[slice_idx]
    neun_slice = neun_stack[slice_idx]
    
    print(f"    TIFF stack shape: {dapi_stack.shape}")
    print(f"    Using slice {slice_idx} (Series {series_num})")
    
    # threshold dapi to find nuclei
    dapi = dapi_slice.astype(float)
    fosb = fosb_slice.astype(float)
    gfap = gfap_slice.astype(float)
    neun = neun_slice.astype(float)
    
    dapi_binary = dapi > DAPI_THRESHOLD
    dapi_binary = clear_border(dapi_binary)
    dapi_binary = remove_small_objects(dapi_binary, min_size=MIN_SIZE)
    
    # load roi exclusion if exists
    roi_points = load_roi_exclusion_mask(replicate_id, 0, roi_dir)
    exclusion_mask = create_exclusion_mask(roi_points, dapi_binary.shape)
    
    if np.any(exclusion_mask):
        dapi_binary = dapi_binary & ~exclusion_mask
    
    # label nuclei
    nuclei_labeled = measure.label(dapi_binary)
    
    # classify as astrocyte or neuron based on gfap
    nuclei_template = []
    astrocyte_count = 0
    
    for region in measure.regionprops(nuclei_labeled):
        coords = region.coords
        
        # get avg gfap intensity
        gfap_vals = [gfap[coord[0], coord[1]] for coord in coords]
        avg_gfap = np.mean(gfap_vals)
        
        is_astrocyte = avg_gfap >= GFAP_THRESHOLD
        
        nuclei_template.append({
            'nucleus_id': region.label,
            'is_astrocyte': is_astrocyte
        })
        
        if is_astrocyte:
            astrocyte_count += 1
    
    total_nuclei = len(nuclei_template)
    print(f"    Template: {total_nuclei} nuclei ({astrocyte_count} astrocytes)")
    
    return {
        'nuclei_labeled': nuclei_labeled,
        'nuclei_template': nuclei_template,
        'total_nuclei': total_nuclei,
        'total_astrocytes': astrocyte_count
    }

def apply_tiff_template_to_nd2_series(replicate_id, nd2_path, tiff_template, target_series_num, output_dir):
    # apply template to nd2 z-slices
    print(f"    Applying TIFF template to ND2 Series {target_series_num}...")
    
    series_slice_data = []
    
    try:
        with ND2Reader(nd2_path) as images:
            images.bundle_axes = ['c', 'y', 'x']
            images.iter_axes = ['z', 'v']
            
            num_z_slices = images.sizes.get('z', 1)
            num_fovs = images.sizes.get('v', 1)
            
            nuclei_labeled = tiff_template['nuclei_labeled']
            nuclei_template = tiff_template['nuclei_template']
            total_nuclei = tiff_template['total_nuclei']
            total_astrocytes = tiff_template['total_astrocytes']
            
            fov_idx = target_series_num - 1
            
            if fov_idx >= num_fovs:
                print(f"    ERROR: Series {target_series_num} doesn't exist")
                return series_slice_data
            
            print(f"    Processing {num_z_slices} z-slices...")
            
            # process all z-slices in this series
            for z_idx in range(num_z_slices):
                frame_idx = fov_idx * num_z_slices + z_idx
                frame = images[frame_idx]
                
                if frame.ndim == 3 and frame.shape[0] == 4:
                    fosb_slice = frame[1]  # fosb is channel 1
                    
                    # analyze this slice
                    slice_data = analyze_nd2_zslice_with_tiff_template(
                        fosb_slice, nuclei_labeled, nuclei_template,
                        replicate_id, target_series_num, z_idx,
                        total_nuclei, total_astrocytes
                    )
                    
                    if slice_data:
                        series_slice_data.append(slice_data)
                    
                    # print progress every 10 slices
                    if z_idx % 10 == 0 and slice_data:
                        print(f"      Z-slice {z_idx}: FosB ratio = {slice_data['fosb_per_pixel_ratio']:.4f}")
            
            print(f"    Done: {len(series_slice_data)} slices analyzed")
    
    except Exception as e:
        print(f"    ERROR: {e}")
    
    return series_slice_data

def analyze_nd2_zslice_with_tiff_template(fosb_slice, nuclei_labeled, nuclei_template,
                                        replicate_id, series_number, z_idx,
                                        total_nuclei, total_astrocytes):
    # measure fosb in each nucleus from template
    all_nuclei_fosb_pixels = []
    astrocyte_fosb_pixels = []
    
    for nucleus_info in nuclei_template:
        nucleus_id = nucleus_info['nucleus_id']
        is_astrocyte = nucleus_info['is_astrocyte']
        
        mask = nuclei_labeled == nucleus_id
        
        if np.any(mask):
            # get fosb pixels for this nucleus
            nucleus_fosb_pixels = fosb_slice[mask]
            all_nuclei_fosb_pixels.extend(nucleus_fosb_pixels.tolist())
            
            if is_astrocyte:
                astrocyte_fosb_pixels.extend(nucleus_fosb_pixels.tolist())
    
    # calculate mean fosb ratio
    mean_fosb_all = np.mean(all_nuclei_fosb_pixels) if all_nuclei_fosb_pixels else 0.0
    mean_fosb_astrocytes = np.mean(astrocyte_fosb_pixels) if astrocyte_fosb_pixels else 0.0
    
    fosb_per_pixel_ratio = (mean_fosb_astrocytes / mean_fosb_all 
                           if mean_fosb_all > 0 else 0.0)
    
    # also calc sum-based metrics
    total_fosb_intensity = np.sum(all_nuclei_fosb_pixels) if all_nuclei_fosb_pixels else 0.0
    astrocyte_fosb_intensity = np.sum(astrocyte_fosb_pixels) if astrocyte_fosb_pixels else 0.0
    
    slice_data = {
        'replicate_id': replicate_id,
        'series_number': series_number,
        'slice_index': z_idx,
        'total_nuclei': total_nuclei,
        'astrocyte_count': total_astrocytes,
        'neuron_count': total_nuclei - total_astrocytes,
        'astrocyte_to_all_nuclei': total_astrocytes / total_nuclei if total_nuclei > 0 else 0,
        'total_fosb_pixels': len(all_nuclei_fosb_pixels),
        'astrocyte_fosb_pixels': len(astrocyte_fosb_pixels),
        'total_fosb_intensity_all_pixels': total_fosb_intensity,
        'astrocyte_fosb_intensity_per_pixel': astrocyte_fosb_intensity,
        'fosb_per_pixel_ratio': fosb_per_pixel_ratio,
        'mean_fosb_per_pixel_all': mean_fosb_all,
        'mean_fosb_per_pixel_astrocytes': mean_fosb_astrocytes,
        'group': 'Control' if '_C' in replicate_id else 'Experimental',
        'data_source': 'ND2_with_TIFF_nuclei_template',
        'template_source': 'TIFF_MIP_exact_fixed_raw_thresholds',
        'nuclei_method': 'EXACT_fixed_raw_threshold_analysis',
        'dapi_threshold': DAPI_THRESHOLD,
        'gfap_threshold': GFAP_THRESHOLD,
        'min_size': MIN_SIZE,
        'analysis_date': datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    }
    
    return slice_data

def main():
    print("TIFF NUCLEI TEMPLATE + ND2 FOSB ANALYSIS")
    print("="*80)
    print(f"DAPI: {DAPI_THRESHOLD}, GFAP: {GFAP_THRESHOLD}, MinSize: {MIN_SIZE}")
    print("="*80)
    
    # setup directories
    tiff_dir = "/Users/charlessander/Desktop/patzke lab computing/Datta Lab/2025-08-31_FOSB_IHC"
    nd2_dir = "/Users/charlessander/Desktop/patzke lab computing/Datta Lab/z stacks/2025-08-31_FOSB_IHC_deconv"
    roi_dir = "/Users/charlessander/Desktop/patzke lab computing/Datta Lab/roi_exclusions"
    output_dir = "/Users/charlessander/Desktop/patzke lab computing/Datta Lab/exact_tiff_nuclei_nd2_fosb_analysis"
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(roi_dir, exist_ok=True)
    
    # create mapping
    print("\nCreating TIFF-to-ND2 mapping...")
    mapping_df = create_tiff_nd2_mapping()
    mapping_path = os.path.join(output_dir, 'tiff_nd2_mapping.csv')
    mapping_df.to_csv(mapping_path, index=False)
    print(f"  Saved mapping: {len(mapping_df)} series")
    
    # load tiff files
    print("\nLoading TIFF data...")
    tiff_replicates = get_biological_replicates_exact_method(tiff_dir)
    print(f"  Found replicates: {list(tiff_replicates.keys())}")
    
    # process each replicate
    print("\nProcessing replicates...")
    all_slice_data = []
    
    nd2_files = [
        ('176130_C', '176130_C.nd2'),
        ('178992_C', '178992_C.nd2'),
        ('178993_C', '178993_C.nd2'),
        ('179041_E', '179041_E.nd2'),
        ('179042_E', '179042_E.nd2')
    ]
    
    for replicate_id, nd2_file in nd2_files:
        tiff_id = replicate_id.split('_')[0]
        nd2_path = os.path.join(nd2_dir, nd2_file)
        
        if not os.path.exists(nd2_path):
            print(f"ERROR: ND2 file not found: {nd2_path}")
            continue
        
        if tiff_id not in tiff_replicates:
            print(f"ERROR: TIFF not found: {tiff_id}")
            continue
        
        group = 'Control' if '_C' in replicate_id else 'Experimental'
        print(f"\n{'='*80}")
        print(f"PROCESSING {replicate_id} ({group})")
        print(f"{'='*80}")
        
        # get series info
        series_info = mapping_df[mapping_df['replicate_id'] == replicate_id]
        num_series = len(series_info)
        print(f"  {num_series} series to process")
        
        slice_data = []
        
        # process each series
        for _, series_row in series_info.iterrows():
            series_num = series_row['series_number']
            print(f"\n  === SERIES {series_num} ===")
            
            # create template for this series
            tiff_template = create_tiff_nuclei_template_exact_method(
                tiff_replicates[tiff_id], replicate_id, series_num, roi_dir, output_dir
            )
            
            if tiff_template is None:
                print(f"  ERROR: Could not create template for Series {series_num}")
                continue
            
            # apply template to nd2 slices
            series_slice_data = apply_tiff_template_to_nd2_series(
                replicate_id, nd2_path, tiff_template, series_num, output_dir
            )
            
            slice_data.extend(series_slice_data)
        
        all_slice_data.extend(slice_data)
    
    # save results
    print(f"\n{'='*80}")
    print("SAVING RESULTS")
    print(f"{'='*80}")
    
    if all_slice_data:
        slice_df = pd.DataFrame(all_slice_data)
        slice_df.to_csv(os.path.join(output_dir, 'exact_tiff_nuclei_nd2_fosb_slice_data.csv'), index=False)
        print(f"  Saved: {len(slice_df)} slices")
        
        # print summary stats
        control_slices = slice_df[slice_df['group'] == 'Control']
        exp_slices = slice_df[slice_df['group'] == 'Experimental']
        
        print(f"\nSUMMARY:")
        print(f"  Total slices: {len(slice_df)}")
        print(f"  Control slices: {len(control_slices)}")
        print(f"  Experimental slices: {len(exp_slices)}")
        
        if len(control_slices) > 0 and len(exp_slices) > 0:
            print(f"  Control mean FosB ratio: {control_slices['fosb_per_pixel_ratio'].mean():.4f}")
            print(f"  Experimental mean FosB ratio: {exp_slices['fosb_per_pixel_ratio'].mean():.4f}")
            
            t_stat, p_value = stats.ttest_ind(exp_slices['fosb_per_pixel_ratio'], control_slices['fosb_per_pixel_ratio'])
            print(f"  t-test p-value: {p_value:.6f}")
            print(f"  Effect size: {exp_slices['fosb_per_pixel_ratio'].mean()/control_slices['fosb_per_pixel_ratio'].mean():.1f}x")
        
        print(f"\nDone! Results saved to: {output_dir}")
    else:
        print("ERROR: No data generated")

if __name__ == "__main__":
    main()
