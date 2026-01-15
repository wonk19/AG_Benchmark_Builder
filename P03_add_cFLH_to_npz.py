"""
P03_add_cFLH_to_npz.py

Add cFLH values to all pixels in archived cluster files.

This script:
1. Reads all cluster_*.npz files in GUI_archive/cluster_archive/
2. Calculates cFLH for every pixel using bio-optical model
3. Saves new files with '_cFLH' suffix containing original data + cFLH values
4. Skips files that already have '_cFLH' version

Note: This process is computationally intensive and may take hours for large datasets.
"""

import sys
import numpy as np
from pathlib import Path
from tqdm import tqdm

# Add paths
base_dir = Path(__file__).parent
sys.path.insert(0, str(base_dir / 'MargModel'))

import MargModel


def load_r2_rank2_parameters():
    """
    Load R2 Rank 2 parameters.
    Uses default values if parameter file not found.
    """
    return {
        'eta': 0.775,
        'g0': 0.00020,
        'nu': 1.000
    }


def calculate_fluorescence(wl, chla):
    """
    Calculate fluorescence spectrum based on Chla.
    """
    if chla <= 0:
        return np.zeros_like(wl)
    
    # Peak wavelength
    p_wl = 4.239291
    q_wl = 678.895348
    peak_wl = p_wl * np.log(chla) + q_wl
    
    # Fluorescence intensity
    a3 = 0.0000102764
    a2 = 0.0000551991
    a1 = 0.0000666462
    a0 = -0.0000137463
    
    log_chla = np.log(chla)
    fl_intensity = a3 * log_chla**3 + a2 * log_chla**2 + a1 * log_chla + a0
    fl_intensity = max(0, fl_intensity)
    
    # Gaussian spectrum
    sigma = 14.862
    fluorescence = fl_intensity * np.exp(-0.5 * ((wl - peak_wl) / sigma)**2)
    
    return fluorescence


def optimize_water_quality(rrs_obs, wl_obs, params, mode='relative'):
    """
    Optimize water quality parameters to match observed Rrs.
    """
    from scipy.optimize import minimize
    
    # Wavelength range for optimization
    wl_fit_range = (wl_obs >= 443) & (wl_obs <= 745)
    wl_fit = wl_obs[wl_fit_range]
    rrs_fit = rrs_obs[wl_fit_range]
    
    # Remove NaN values
    valid = np.isfinite(rrs_fit)
    wl_fit = wl_fit[valid]
    rrs_fit = rrs_fit[valid]
    
    if len(wl_fit) < 2:
        return {'success': False, 'message': 'Insufficient valid data points'}
    
    def objective(x):
        """Objective function based on selected mode"""
        chla, mineral, acdom = x
        
        # Bounds check
        if chla < 0.1 or chla > 1000:
            return 1e10
        if mineral < 0 or mineral > 50:
            return 1e10
        if acdom < 0 or acdom > 5:
            return 1e10
        
        try:
            # Run forward model
            rrs_model = MargModel.runForward(wl_fit, chla, mineral, acdom)
            
            # Add fluorescence
            fluorescence = calculate_fluorescence(wl_fit, chla)
            rrs_model_with_fluo = rrs_model + fluorescence
            
            if mode == 'strict':
                # STRICT MODE: Match absolute Rrs values directly
                errors = []
                for i, wl in enumerate(wl_fit):
                    error = (rrs_model_with_fluo[i] - rrs_fit[i])**2
                    
                    if wl == 709:
                        errors.append(error * 10.0)
                    elif wl == 745:
                        errors.append(error * 5.0)
                    else:
                        errors.append(error * 1.0)
                
                rmse = np.sqrt(np.mean(errors))
                return rmse
                
            else:  # mode == 'relative'
                # RELATIVE MODE: Match spectral differences
                errors = []
                
                # Match absolute values at 709, 745
                for i, wl in enumerate(wl_fit):
                    if wl == 709:
                        error = (rrs_model_with_fluo[i] - rrs_fit[i])**2
                        errors.append(error * 10.0)
                    elif wl == 745:
                        error = (rrs_model_with_fluo[i] - rrs_fit[i])**2
                        errors.append(error * 5.0)
                
                # Match spectral differences for other bands
                for i in range(1, len(wl_fit)):
                    obs_diff = rrs_fit[i] - rrs_fit[i-1]
                    model_diff = rrs_model_with_fluo[i] - rrs_model_with_fluo[i-1]
                    error = (model_diff - obs_diff)**2
                    
                    if wl_fit[i-1] == 680 and wl_fit[i] == 709:
                        errors.append(error * 3.0)
                    else:
                        errors.append(error * 1.0)
                
                if len(errors) == 0:
                    return 1e10
                
                rmse = np.sqrt(np.mean(errors))
                return rmse
            
        except:
            return 1e10
    
    # Bounds
    bounds = [(0.1, 1000), (0, 50), (0, 5)]
    
    # Try multiple initial guesses
    initial_guesses = []
    for chla in [30, 50, 70, 100, 150, 200]:
        for mineral in [0.1, 0.5, 1.0, 3.0]:
            for acdom in [0.1, 0.3, 0.7]:
                initial_guesses.append([chla, mineral, acdom])
    
    best_result = None
    best_rmse = 1e10
    
    for x0 in initial_guesses:
        try:
            result = minimize(objective, x0, method='L-BFGS-B', bounds=bounds,
                            options={'maxiter': 300, 'ftol': 1e-9})
            
            if result.success and result.fun < best_rmse:
                best_rmse = result.fun
                best_result = result
        except:
            continue
        
        if best_rmse < 0.0005:
            break
    
    if best_result is not None and best_result.success:
        chla_opt, mineral_opt, acdom_opt = best_result.x
        rmse = best_result.fun
        
        return {
            'success': True,
            'chla': chla_opt,
            'mineral': mineral_opt,
            'acdom': acdom_opt,
            'rmse': rmse
        }
    else:
        return {'success': False, 'message': 'All optimization attempts failed'}


def calculate_cflh_for_pixel(wl_obs, rrs_obs, params, mode='relative'):
    """
    Calculate cFLH for a single pixel.
    
    Returns:
    --------
    cflh_value : float
        cFLH value, or np.nan if calculation failed
    """
    try:
        # Optimize water quality parameters
        opt_result = optimize_water_quality(rrs_obs, wl_obs, params, mode=mode)
        
        if not opt_result['success']:
            return np.nan
        
        chla = opt_result['chla']
        mineral = opt_result['mineral']
        acdom = opt_result['acdom']
        
        # Generate hyperspectral Rrs with mineral=0
        wl_hyper = np.arange(400., 749.5, 0.5)
        
        # Without mineral (corrected)
        rrs_no_mineral = MargModel.runForward(wl_hyper, chla, 0.0, acdom)
        fluo_no_mineral = calculate_fluorescence(wl_hyper, chla)
        rrs_total_no_mineral = rrs_no_mineral + fluo_no_mineral
        
        # Find peak between 660-745nm
        idx_660 = np.argmin(np.abs(wl_hyper - 660))
        idx_745 = np.argmin(np.abs(wl_hyper - 745))
        
        idx_range = (wl_hyper >= 660) & (wl_hyper <= 745)
        wl_range = wl_hyper[idx_range]
        rrs_range = rrs_total_no_mineral[idx_range]
        
        # Baseline interpolation
        baseline_range = rrs_total_no_mineral[idx_660] + \
                        (rrs_total_no_mineral[idx_745] - rrs_total_no_mineral[idx_660]) * \
                        (wl_range - 660) / (745 - 660)
        
        # Height above baseline
        height_above_baseline = rrs_range - baseline_range
        idx_peak_local = np.argmax(height_above_baseline)
        cFLH = height_above_baseline[idx_peak_local]
        
        return cFLH
        
    except Exception as e:
        return np.nan


def process_cluster_file(input_file, output_file, mode='relative'):
    """
    Process a single cluster file and add cFLH values.
    """
    print(f"\nProcessing: {input_file.name}")
    
    # Load input file
    data = np.load(input_file, allow_pickle=True)
    clusters = data['archived_clusters'].item()
    
    print(f"  Found {len(clusters)} clusters")
    
    # Load parameters once
    params = load_r2_rank2_parameters()
    
    # GOCI-II wavelengths
    wavelengths = np.array([380, 412, 443, 490, 510, 555, 620, 660, 680, 709, 745, 865])
    
    # Process each cluster
    total_pixels = sum(cluster_data['n_pixels'] for cluster_data in clusters.values())
    
    print(f"  Total pixels to process: {total_pixels}")
    print(f"  Mode: {mode}")
    
    with tqdm(total=total_pixels, desc="  Calculating cFLH", ncols=80) as pbar:
        for cluster_id, cluster_data in clusters.items():
            n_pixels = cluster_data['n_pixels']
            
            # Initialize cFLH array for this cluster
            cflh_values = np.zeros(n_pixels)
            
            # Process each pixel
            for pixel_idx in range(n_pixels):
                # Extract Rrs for this pixel
                rrs_pixel = np.array([
                    cluster_data['Rrs_380'][pixel_idx],
                    cluster_data['Rrs_412'][pixel_idx],
                    cluster_data['Rrs_443'][pixel_idx],
                    cluster_data['Rrs_490'][pixel_idx],
                    cluster_data['Rrs_510'][pixel_idx],
                    cluster_data['Rrs_555'][pixel_idx],
                    cluster_data['Rrs_620'][pixel_idx],
                    cluster_data['Rrs_660'][pixel_idx],
                    cluster_data['Rrs_680'][pixel_idx],
                    cluster_data['Rrs_709'][pixel_idx],
                    cluster_data['Rrs_745'][pixel_idx],
                    cluster_data['Rrs_865'][pixel_idx]
                ])
                
                # Calculate cFLH
                cflh_value = calculate_cflh_for_pixel(wavelengths, rrs_pixel, params, mode=mode)
                cflh_values[pixel_idx] = cflh_value
                
                pbar.update(1)
            
            # Add cFLH to cluster data
            cluster_data['cFLH'] = cflh_values
    
    # Save new file
    np.savez_compressed(output_file, archived_clusters=clusters)
    print(f"  Saved: {output_file.name}")
    
    # Print statistics
    all_cflh = []
    all_flh = []
    for cluster_data in clusters.values():
        all_cflh.extend(cluster_data['cFLH'][np.isfinite(cluster_data['cFLH'])])
        all_flh.extend(cluster_data['FLH'])
    
    if len(all_cflh) > 0:
        print(f"  cFLH statistics:")
        print(f"    Valid pixels: {len(all_cflh)} / {total_pixels} ({len(all_cflh)/total_pixels*100:.1f}%)")
        print(f"    Range: [{np.min(all_cflh):.6f}, {np.max(all_cflh):.6f}]")
        print(f"    Mean: {np.mean(all_cflh):.6f}")
        print(f"  FLH statistics:")
        print(f"    Range: [{np.min(all_flh):.6f}, {np.max(all_flh):.6f}]")
        print(f"    Mean: {np.mean(all_flh):.6f}")
        print(f"  cFLH/FLH ratio: {np.mean(all_cflh)/np.mean(all_flh):.3f}")


def main():
    print("="*70)
    print("P03: Add cFLH to Archived Cluster Files")
    print("="*70)
    
    # Select mode
    MODE = 'relative'  # Change to 'strict' if needed
    print(f"\nOptimization mode: {MODE}")
    
    # Archive directory
    base_dir = Path(__file__).parent
    archive_dir = base_dir / 'results' / 'GUI_archive' / 'cluster_archive'
    
    if not archive_dir.exists():
        print(f"\nError: Archive directory not found: {archive_dir}")
        return
    
    # Find all cluster files
    cluster_files = sorted(archive_dir.glob("clusters_*.npz"))
    
    # Filter out files that already have '_cFLH' suffix
    files_to_process = []
    for f in cluster_files:
        if '_cFLH' in f.stem:
            print(f"Skipping (already processed): {f.name}")
            continue
        
        # Check if output file already exists
        output_file = f.parent / f"{f.stem}_cFLH.npz"
        if output_file.exists():
            print(f"Skipping (output exists): {f.name} -> {output_file.name}")
            continue
        
        files_to_process.append(f)
    
    if not files_to_process:
        print("\nNo files to process!")
        return
    
    print(f"\nFound {len(files_to_process)} files to process:")
    for f in files_to_process:
        print(f"  - {f.name}")
    
    # Process each file
    print("\n" + "="*70)
    print("Starting processing...")
    print("="*70)
    
    for i, input_file in enumerate(files_to_process):
        print(f"\n[{i+1}/{len(files_to_process)}]")
        
        output_file = input_file.parent / f"{input_file.stem}_cFLH.npz"
        
        try:
            process_cluster_file(input_file, output_file, mode=MODE)
        except Exception as e:
            print(f"  ERROR processing {input_file.name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n" + "="*70)
    print("Processing complete!")
    print("="*70)


if __name__ == "__main__":
    main()

