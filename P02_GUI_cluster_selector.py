"""
P02_GUI_cluster_selector.py

GOCI-II Clustering GUI
- Load cropped GOCI-II data
- Perform online clustering (Bayesian GMM + spatial splitting)
- View and interact with clusters
- Archive selected clusters
- Manage global cluster database
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from pathlib import Path
from datetime import datetime
import tkinter as tk
from tkinter import ttk, messagebox
import json
import pickle
import pandas as pd
from scipy.linalg import sqrtm
from scipy.ndimage import label as ndimage_label, sobel
from sklearn.mixture import BayesianGaussianMixture

# Add library path
sys.path.insert(0, str(Path(__file__).parent / 'library'))

from goci2_reader import get_flag_by_name


def compute_ndvi(band_865, band_660):
    return (band_865 - band_660) / (band_865 + band_660 + 1e-8)


def compute_flh(band_660, band_709, band_745):
    """Compute FLH with baseline 660-745nm"""
    lambda_660, lambda_709, lambda_745 = 660.0, 709.0, 745.0
    baseline = band_660 + (lambda_709 - lambda_660) / (lambda_745 - lambda_660) * (band_745 - band_660)
    return band_709 - baseline


def compute_mnf_scene(X, n_components=4):
    cov_data = np.cov(X.T)
    diffs = np.diff(X, axis=0)
    cov_noise = np.cov(diffs.T) / 2.0
    cov_noise_reg = cov_noise + 1e-6 * np.eye(cov_noise.shape[0])
    
    try:
        L = np.linalg.cholesky(cov_noise_reg)
    except np.linalg.LinAlgError:
        eigenvalues_noise, eigenvectors_noise = np.linalg.eigh(cov_noise_reg)
        eigenvalues_noise = np.maximum(eigenvalues_noise, 1e-10)
        L = eigenvectors_noise @ np.diag(np.sqrt(eigenvalues_noise)) @ eigenvectors_noise.T
    
    L_inv = np.linalg.inv(L)
    cov_data_whitened = L_inv @ cov_data @ L_inv.T
    eigenvalues, eigenvectors = np.linalg.eigh(cov_data_whitened)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    transformation = L_inv.T @ eigenvectors
    Z = X @ transformation[:, :n_components]
    
    return Z, transformation[:, :n_components], eigenvalues[:n_components]


def robust_scale(F):
    median = np.median(F, axis=0)
    mad = np.median(np.abs(F - median), axis=0)
    mad = np.where(mad < 1e-10, 1.0, mad)
    return (F - median) / (1.4826 * mad)


def l2_normalize(F):
    norms = np.linalg.norm(F, axis=1, keepdims=True)
    norms = np.where(norms < 1e-10, 1.0, norms)
    return F / norms


def weighted_mean(X, w):
    w = w / (w.sum() + 1e-12)
    return (X.T @ w).ravel()


def weighted_cov(X, w):
    w = w / (w.sum() + 1e-12)
    mu = weighted_mean(X, w)
    X_centered = X - mu
    return (X_centered.T * w) @ X_centered


def shrink_cov(Sigma, lam=0.1):
    d = Sigma.shape[0]
    I = np.eye(d)
    trace_sigma = np.trace(Sigma)
    return (1 - lam) * Sigma + lam * (trace_sigma / d) * I


def spectral_angle(mu1, mu2):
    norm1 = np.linalg.norm(mu1)
    norm2 = np.linalg.norm(mu2)
    if norm1 < 1e-10 or norm2 < 1e-10:
        return np.pi / 2
    mu1_normalized = mu1 / norm1
    mu2_normalized = mu2 / norm2
    cos_angle = np.dot(mu1_normalized, mu2_normalized)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    return np.arccos(cos_angle)


def symKL_gaussian(mu1, Sigma1, mu2, Sigma2):
    d = len(mu1)
    try:
        Sigma1_inv = np.linalg.inv(Sigma1 + 1e-6 * np.eye(d))
        Sigma2_inv = np.linalg.inv(Sigma2 + 1e-6 * np.eye(d))
    except:
        return 1e10
    diff = mu1 - mu2
    kl12 = 0.5 * (np.trace(Sigma2_inv @ Sigma1) + diff.T @ Sigma2_inv @ diff - d + np.log(np.linalg.det(Sigma2) / np.linalg.det(Sigma1)))
    kl21 = 0.5 * (np.trace(Sigma1_inv @ Sigma2) + diff.T @ Sigma1_inv @ diff - d + np.log(np.linalg.det(Sigma1) / np.linalg.det(Sigma2)))
    return kl12 + kl21


def split_clusters_by_spatial_connectivity(labels_flat, water_mask, H, W, min_size=10):
    labels_2d = np.full((H, W), -1, dtype=np.int32)
    labels_2d_flat = labels_2d.ravel()
    labels_2d_flat[water_mask] = labels_flat
    labels_2d = labels_2d_flat.reshape(H, W)
    
    unique_labels = np.unique(labels_flat)
    unique_labels = unique_labels[unique_labels >= 0]
    
    new_label_counter = 0
    new_labels_2d = np.full((H, W), -1, dtype=np.int32)
    spatial_to_spectral_map = {}
    
    for original_label in unique_labels:
        mask = (labels_2d == original_label)
        labeled_components, num_components = ndimage_label(mask, structure=np.ones((3, 3)))
        
        largest_component_label = None
        largest_component_size = 0
        
        for component_id in range(1, num_components + 1):
            component_mask = (labeled_components == component_id)
            component_size = component_mask.sum()
            
            if component_size >= min_size:
                new_labels_2d[component_mask] = new_label_counter
                spatial_to_spectral_map[new_label_counter] = original_label
                if component_size > largest_component_size:
                    largest_component_size = component_size
                    largest_component_label = new_label_counter
                new_label_counter += 1
        
        for component_id in range(1, num_components + 1):
            component_mask = (labeled_components == component_id)
            component_size = component_mask.sum()
            
            if component_size < min_size:
                if largest_component_label is not None:
                    new_labels_2d[component_mask] = largest_component_label
                else:
                    new_labels_2d[component_mask] = new_label_counter
                    spatial_to_spectral_map[new_label_counter] = original_label
                    new_label_counter += 1
    
    new_labels_2d_flat = new_labels_2d.ravel()
    new_labels_flat = new_labels_2d_flat[water_mask]
    
    return new_labels_flat, new_label_counter, spatial_to_spectral_map


class GlobalClusterManager:
    def __init__(self, n_bands=12):
        self.clusters = []
        self.n_bands = n_bands
        self.cluster_id_counter = 0
        
    def add_new_cluster(self, X, w, meta=None):
        mu = weighted_mean(X, w)
        mu_norm = np.linalg.norm(mu)
        if mu_norm > 1e-10:
            mu = mu / mu_norm
        Sigma = weighted_cov(X, w) + 1e-5 * np.eye(self.n_bands)
        Sigma = shrink_cov(Sigma, lam=0.1)
        cluster = {
            'id': self.cluster_id_counter,
            'mu': mu,
            'Sigma': Sigma,
            'n_samples': len(X),
            'weight_sum': w.sum(),
            'meta': meta or {}
        }
        self.clusters.append(cluster)
        self.cluster_id_counter += 1
        return cluster['id']
    
    def update_cluster(self, k, X, w, rho=0.1):
        if k >= len(self.clusters):
            return
        cluster = self.clusters[k]
        mu_new = weighted_mean(X, w)
        mu_new_norm = np.linalg.norm(mu_new)
        if mu_new_norm > 1e-10:
            mu_new = mu_new / mu_new_norm
        Sigma_new = weighted_cov(X, w) + 1e-5 * np.eye(self.n_bands)
        Sigma_new = shrink_cov(Sigma_new, lam=0.1)
        cluster['mu'] = (1 - rho) * cluster['mu'] + rho * mu_new
        cluster['mu'] = cluster['mu'] / (np.linalg.norm(cluster['mu']) + 1e-10)
        cluster['Sigma'] = (1 - rho) * cluster['Sigma'] + rho * Sigma_new
        cluster['n_samples'] += len(X)
        cluster['weight_sum'] += w.sum()
    
    def find_best_match(self, X, w, sam_threshold=0.15, kl_threshold=5.0):
        if len(self.clusters) == 0:
            return None, None
        mu_local = weighted_mean(X, w)
        mu_local_norm = np.linalg.norm(mu_local)
        if mu_local_norm > 1e-10:
            mu_local = mu_local / mu_local_norm
        Sigma_local = weighted_cov(X, w) + 1e-5 * np.eye(self.n_bands)
        Sigma_local = shrink_cov(Sigma_local, lam=0.1)
        
        best_k = None
        best_score = np.inf
        
        for k, cluster in enumerate(self.clusters):
            sam = spectral_angle(mu_local, cluster['mu'])
            kl = symKL_gaussian(mu_local, Sigma_local, cluster['mu'], cluster['Sigma'])
            
            if sam < sam_threshold and kl < kl_threshold:
                score = sam + 0.1 * kl
                if score < best_score:
                    best_score = score
                    best_k = k
        
        return best_k, best_score


class ClusteringGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("GOCI-II Clustering Viewer")
        self.root.geometry("1800x900")
        
        # Paths
        self.base_dir = Path(__file__).parent
        self.cropped_dir = self.base_dir / 'results' / 'GUI_archive' / 'cropped_data'
        self.nfis_dir = self.base_dir / 'data' / 'nfis'
        self.archive_dir = self.base_dir / 'results' / 'GUI_archive'
        self.clustering_results_dir = self.archive_dir / 'clustering_results'
        self.clustering_results_dir.mkdir(parents=True, exist_ok=True)
        self.cluster_archive_dir = self.archive_dir / 'cluster_archive'
        self.cluster_archive_dir.mkdir(parents=True, exist_ok=True)
        
        self.current_data = None
        self.current_file = None
        self.nfis_points = None
        
        self.global_manager = None
        self.scene_results = []
        self.current_labels = None
        self.current_global_labels = None
        self.current_flh = None
        self.current_flh_2d = None
        self.current_global_flh_map = None
        self.current_Rrs_all = None
        self.water_mask = None
        
        self.show_flh = False
        self.hold_spectra = False
        self.spectra_list = []
        self.highlighted_clusters = set()
        self.delete_mode = False
        
        self.wavelengths = [380, 412, 443, 490, 510, 555, 620, 660, 680, 709, 745, 865]
        
        self.load_global_clusters()
        self.create_widgets()
        self.load_dates()
        
    def load_global_clusters(self):
        global_clusters_file = self.archive_dir / 'global_clusters.pkl'
        scene_results_file = self.archive_dir / 'scene_results.json'
        
        if global_clusters_file.exists():
            with open(global_clusters_file, 'rb') as f:
                loaded_data = pickle.load(f)
            
            if isinstance(loaded_data, GlobalClusterManager):
                self.global_manager = loaded_data
            elif isinstance(loaded_data, list):
                self.global_manager = GlobalClusterManager(n_bands=12)
                self.global_manager.clusters = loaded_data
                if loaded_data:
                    max_id = max(c['id'] for c in loaded_data)
                    self.global_manager.cluster_id_counter = max_id + 1
            else:
                self.global_manager = GlobalClusterManager(n_bands=12)
            
            print(f"Loaded {len(self.global_manager.clusters)} global clusters")
        else:
            self.global_manager = GlobalClusterManager(n_bands=12)
            print("Initialized new global cluster manager")
        
        if scene_results_file.exists():
            with open(scene_results_file, 'r') as f:
                self.scene_results = json.load(f)
            print(f"Loaded {len(self.scene_results)} scene results")
        
    def save_global_clusters(self):
        global_clusters_file = self.archive_dir / 'global_clusters.pkl'
        scene_results_file = self.archive_dir / 'scene_results.json'
        
        with open(global_clusters_file, 'wb') as f:
            pickle.dump(self.global_manager, f)
        
        with open(scene_results_file, 'w') as f:
            json.dump(self.scene_results, f, indent=2)
        
        print("Saved global clusters and scene results")
    
    def create_widgets(self):
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.TOP, fill=tk.X)
        
        ttk.Label(control_frame, text="Date:").pack(side=tk.LEFT, padx=(0, 5))
        self.date_combo = ttk.Combobox(control_frame, state='readonly', width=15)
        self.date_combo.pack(side=tk.LEFT, padx=(0, 10))
        self.date_combo.bind('<<ComboboxSelected>>', self.on_date_change)
        
        ttk.Label(control_frame, text="File:").pack(side=tk.LEFT, padx=(0, 5))
        self.file_combo = ttk.Combobox(control_frame, state='readonly', width=60)
        self.file_combo.pack(side=tk.LEFT, padx=(0, 10))
        self.file_combo.bind('<<ComboboxSelected>>', self.on_file_select)
        
        self.cluster_button = ttk.Button(control_frame, text="Clustering", 
                                        command=self.perform_clustering, state=tk.DISABLED)
        self.cluster_button.pack(side=tk.LEFT, padx=5)
        
        self.hold_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(control_frame, text="Hold", variable=self.hold_var).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(control_frame, text="Clear", command=self.on_clear_spectra).pack(side=tk.LEFT, padx=5)
        
        self.archive_button = ttk.Button(control_frame, text="Archive", 
                                         command=self.archive_selected_clusters, state=tk.DISABLED)
        self.archive_button.pack(side=tk.LEFT, padx=5)
        
        self.display_archived_button = ttk.Button(control_frame, text="Display Archived", 
                                                  command=self.display_archived_clusters, state=tk.DISABLED)
        self.display_archived_button.pack(side=tk.LEFT, padx=5)
        
        self.delete_cluster_button = ttk.Button(control_frame, text="Delete Cluster", 
                                                command=self.toggle_delete_mode, state=tk.DISABLED)
        self.delete_cluster_button.pack(side=tk.LEFT, padx=5)
        
        self.info_label = ttk.Label(control_frame, text="No file loaded")
        self.info_label.pack(side=tk.LEFT, padx=10)
        
        fig_frame = ttk.Frame(self.root)
        fig_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.fig = Figure(figsize=(18, 9))
        self.ax_left = self.fig.add_axes([0.02, 0.05, 0.58, 0.9])
        self.ax_right = self.fig.add_axes([0.65, 0.1, 0.32, 0.8])
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=fig_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        self.canvas.mpl_connect('button_press_event', self.on_click)
    
    def load_dates(self):
        files = sorted(self.cropped_dir.glob('*.npz'))
        dates = set()
        for f in files:
            parts = f.stem.split('_')
            if len(parts) >= 4:
                date = parts[3]
                dates.add(date)
        
        self.dates = sorted(dates)
        self.date_combo['values'] = self.dates
        if self.dates:
            self.date_combo.set(self.dates[0])
            self.on_date_change(None)
    
    def on_date_change(self, event):
        selected_date = self.date_combo.get()
        files = sorted(self.cropped_dir.glob(f'*_{selected_date}_*.npz'))
        
        self.current_files = files
        file_names = [f.name for f in files]
        self.file_combo['values'] = file_names
        
        if file_names:
            self.file_combo.set(file_names[0])
            self.on_file_select(None)
        
        self.load_nfis_for_date(selected_date)
    
    def load_nfis_for_date(self, date):
        nfis_file = self.nfis_dir / f'red_tide_locations_{date[:4]}-{date[4:6]}-{date[6:]}_unique_coords.csv'
        
        if nfis_file.exists():
            df = pd.read_csv(nfis_file)
            self.nfis_points = df[['latitude', 'longitude']].values
            print(f"Loaded {len(self.nfis_points)} NFIS points for {date}")
        else:
            self.nfis_points = None
            print(f"No NFIS file found for {date}")
    
    def on_file_select(self, event):
        selected_file = self.file_combo.get()
        if not selected_file:
            return
        
        filepath = self.cropped_dir / selected_file
        self.load_data(filepath)
    
    def load_data(self, filepath):
        try:
            self.current_data = np.load(filepath)
            self.current_file = filepath
            self.current_labels = None
            self.current_global_labels = None
            self.show_flh = False
            self.spectra_list = []
            self.highlighted_clusters = set()
            
            self.cluster_button.config(state=tk.NORMAL)
            self.display_rrs()
            
            info = f"Loaded: {filepath.name}"
            self.info_label.config(text=info)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")
    
    def display_rrs(self):
        if self.current_data is None:
            return
        
        r_band = self.current_data['Rrs_660']
        g_band = self.current_data['Rrs_555']
        b_band = self.current_data['Rrs_443']
        
        r_norm = self.normalize_band(r_band, 0.0, 0.02)
        g_norm = self.normalize_band(g_band, 0.0, 0.02)
        b_norm = self.normalize_band(b_band, 0.0, 0.02)
        
        rgb = np.dstack([r_norm, g_norm, b_norm])
        
        self.fig.clear()
        self.ax_left = self.fig.add_axes([0.02, 0.05, 0.58, 0.9])
        self.ax_right = self.fig.add_axes([0.65, 0.1, 0.32, 0.8])
        
        self.ax_left.imshow(rgb)
        self.ax_left.axis('off')
        self.ax_left.set_title(f'Rrs RGB - {self.current_file.name}', fontsize=10)
        
        self.ax_right.text(0.5, 0.5, 'Click Clustering button', 
                          ha='center', va='center', transform=self.ax_right.transAxes, fontsize=10)
        self.ax_right.axis('off')
        
        self.plot_nfis_points()
        self.canvas.draw()
    
    def normalize_band(self, band, vmin, vmax):
        valid_mask = np.isfinite(band)
        normalized = np.clip((band - vmin) / (vmax - vmin), 0, 1)
        normalized[~valid_mask] = 0
        return normalized
    
    def plot_nfis_points(self):
        if self.nfis_points is None or self.current_data is None:
            return
        
        lat_grid = self.current_data['latitude']
        lon_grid = self.current_data['longitude']
        H, W = lat_grid.shape
        
        lat_min, lat_max = lat_grid.min(), lat_grid.max()
        lon_min, lon_max = lon_grid.min(), lon_grid.max()
        
        pixel_rows = []
        pixel_cols = []
        
        for lat, lon in self.nfis_points:
            if not (lat_min <= lat <= lat_max and lon_min <= lon <= lon_max):
                continue
            
            distances = np.sqrt((lat_grid - lat)**2 + (lon_grid - lon)**2)
            min_idx = np.unravel_index(np.argmin(distances), distances.shape)
            row, col = min_idx[0], min_idx[1]
            
            if 0 <= row < H and 0 <= col < W:
                pixel_rows.append(row)
                pixel_cols.append(col)
        
        if pixel_rows:
            self.ax_left.scatter(pixel_cols, pixel_rows,
                                c='red', s=30, marker='o',
                                alpha=0.7, edgecolors='darkred', linewidths=1,
                                zorder=10)
    
    def perform_clustering(self):
        if self.current_data is None:
            return
        
        try:
            result_file = self.clustering_results_dir / f'{self.current_file.stem}_labels.npz'
            
            if result_file.exists():
                self.info_label.config(text="Loading existing results...")
                self.root.update()
                self.load_existing_results(result_file)
            else:
                self.info_label.config(text="Clustering in progress...")
                self.root.update()
                self.run_new_clustering(result_file)
            
            self.update_display()
            
            date_str = self.date_combo.get()
            archive_file = self.cluster_archive_dir / f'clusters_{date_str}.npz'
            if archive_file.exists():
                self.display_archived_button.config(state=tk.NORMAL)
                self.delete_cluster_button.config(state=tk.NORMAL)
            
        except Exception as e:
            messagebox.showerror("Error", f"Clustering failed:\n{str(e)}")
            import traceback
            traceback.print_exc()
    
    def load_existing_results(self, result_file):
        result_data = np.load(result_file)
        labels = result_data['labels']
        
        H, W = int(result_data['H']), int(result_data['W'])
        
        current_H, current_W = self.current_data['Rrs_660'].shape
        if H != current_H or W != current_W:
            print(f"Warning: Saved result size ({H}x{W}) != current data size ({current_H}x{current_W})")
            print("Re-running clustering...")
            self.run_new_clustering(result_file)
            return
        
        if 'global_labels' in result_data:
            global_labels = result_data['global_labels']
        else:
            scene_idx = self.find_scene_index(self.current_file.name)
            if scene_idx is not None and scene_idx < len(self.scene_results):
                scene_result = self.scene_results[scene_idx]
                local_to_global = {int(k): int(v) for k, v in scene_result.get('local_to_global', {}).items()}
                
                global_labels_1d = np.full(len(labels), -1, dtype=np.int32)
                for local_id, global_id in local_to_global.items():
                    mask = (labels == local_id)
                    global_labels_1d[mask] = global_id
                
                global_labels = global_labels_1d
            else:
                global_labels = labels
        
        self.water_mask = result_data['water_mask']
        
        self.current_labels = labels.reshape(H, W)
        self.current_global_labels = global_labels.reshape(H, W)
        
        self.prepare_visualization_data(H, W)
        
        n_local = len(np.unique(self.current_labels[self.current_labels >= 0]))
        n_global = len(np.unique(self.current_global_labels[self.current_global_labels >= 0]))
        info = f"Loaded: {n_local} local, {n_global} global clusters"
        self.info_label.config(text=info)
    
    def run_new_clustering(self, result_file):
        labels_2d, global_labels_2d, water_mask, H, W = self.online_clustering()
        
        self.current_labels = labels_2d
        self.current_global_labels = global_labels_2d
        self.water_mask = water_mask
        
        np.savez_compressed(
            result_file,
            labels=labels_2d.ravel(),
            global_labels=global_labels_2d.ravel(),
            water_mask=water_mask,
            H=H,
            W=W
        )
        
        self.prepare_visualization_data(H, W)
        
        n_local = len(np.unique(self.current_labels[self.current_labels >= 0]))
        n_global = len(np.unique(self.current_global_labels[self.current_global_labels >= 0]))
        
        info = f"Complete: {n_local} local, {n_global} global clusters"
        self.info_label.config(text=info)
        
        self.save_global_clusters()
    
    def prepare_visualization_data(self, H, W):
        self.current_Rrs_all = np.stack([self.current_data[f'Rrs_{wl}'].ravel()[self.water_mask] 
                                        for wl in self.wavelengths], axis=-1)
        
        band_660 = self.current_data['Rrs_660'].ravel()[self.water_mask]
        band_709 = self.current_data['Rrs_709'].ravel()[self.water_mask]
        band_745 = self.current_data['Rrs_745'].ravel()[self.water_mask]
        
        self.current_flh = compute_flh(band_660, band_709, band_745)
        
        self.current_flh_2d = np.full((H, W), np.nan)
        flh_2d_flat = self.current_flh_2d.ravel()
        flh_2d_flat[self.water_mask] = self.current_flh
        self.current_flh_2d = flh_2d_flat.reshape(H, W)
        
        global_labels_flat = self.current_global_labels.ravel()[self.water_mask]
        unique_global = np.unique(global_labels_flat[global_labels_flat >= 0])
        
        global_cluster_flh = {}
        for g_id in unique_global:
            mask_g = (global_labels_flat == g_id)
            if mask_g.sum() > 0:
                global_cluster_flh[g_id] = self.current_flh[mask_g].mean()
        
        global_flh_map_flat = np.full(H * W, np.nan)
        for g_id, mean_flh_val in global_cluster_flh.items():
            mask_g = (global_labels_flat == g_id)
            water_indices = np.where(self.water_mask)[0]
            global_indices = water_indices[mask_g]
            global_flh_map_flat[global_indices] = mean_flh_val
        
        self.current_global_flh_map = global_flh_map_flat.reshape(H, W)
    
    def online_clustering(self):
        H, W = self.current_data['Rrs_660'].shape
        
        X_bands = []
        for wl in self.wavelengths:
            band = self.current_data[f'Rrs_{wl}'].ravel()
            X_bands.append(band)
        X = np.column_stack(X_bands)
        
        flag = self.current_data['flag']
        land_mask = get_flag_by_name(flag, 'LAND').ravel().astype(bool)
        cloud_mask = get_flag_by_name(flag, 'CLOUD').ravel().astype(bool)
        
        valid_mask = np.all(np.isfinite(X), axis=1)
        range_mask = (X > -0.01).all(axis=1) & (X < 0.15).all(axis=1)
        
        water_mask = ~land_mask & ~cloud_mask & valid_mask & range_mask
        
        X_water = X[water_mask]
        
        band_660 = X_water[:, self.wavelengths.index(660)]
        band_680 = X_water[:, self.wavelengths.index(680)]
        band_709 = X_water[:, self.wavelengths.index(709)]
        band_745 = X_water[:, self.wavelengths.index(745)]
        band_865 = X_water[:, self.wavelengths.index(865)]
        
        ndvi = compute_ndvi(band_865, band_660)
        flh = compute_flh(band_660, band_709, band_745)
        
        F_spatial = np.column_stack([ndvi, flh])
        Z_mnf, _, _ = compute_mnf_scene(X_water, n_components=4)
        
        F_combined = np.column_stack([Z_mnf, F_spatial])
        F_scaled = robust_scale(F_combined)
        F_norm = l2_normalize(F_scaled)
        
        K_local = 50
        bgmm = BayesianGaussianMixture(
            n_components=K_local,
            covariance_type='full',
            weight_concentration_prior_type='dirichlet_process',
            weight_concentration_prior=0.1,
            reg_covar=1e-6,
            init_params='kmeans',
            max_iter=200,
            tol=1e-3,
            n_init=1,
            random_state=42
        )
        
        bgmm.fit(F_norm)
        labels_spectral = bgmm.predict(F_norm)
        
        labels_spatial, n_spatial, spatial_to_spectral = split_clusters_by_spatial_connectivity(
            labels_spectral, water_mask, H, W, min_size=10
        )
        
        global_labels = np.full(len(labels_spatial), -1, dtype=np.int32)
        local_to_global = {}
        
        for local_id in range(n_spatial):
            mask = (labels_spatial == local_id)
            if mask.sum() == 0:
                continue
            
            X_local = X_water[mask]
            w_local = np.ones(len(X_local))
            
            best_k, best_score = self.global_manager.find_best_match(
                X_local, w_local, sam_threshold=0.15, kl_threshold=5.0
            )
            
            if best_k is not None:
                self.global_manager.update_cluster(best_k, X_local, w_local, rho=0.1)
                global_id = self.global_manager.clusters[best_k]['id']
            else:
                global_id = self.global_manager.add_new_cluster(X_local, w_local)
            
            global_labels[mask] = global_id
            local_to_global[local_id] = global_id
        
        scene_result = {
            'filename': self.current_file.name,
            'n_local_clusters': n_spatial,
            'n_global_clusters': len(self.global_manager.clusters),
            'local_to_global': local_to_global
        }
        self.scene_results.append(scene_result)
        
        labels_2d = np.full((H, W), -1, dtype=np.int32)
        labels_2d_flat = labels_2d.ravel()
        labels_2d_flat[water_mask] = labels_spatial
        labels_2d = labels_2d_flat.reshape(H, W)
        
        global_labels_2d = np.full((H, W), -1, dtype=np.int32)
        global_labels_2d_flat = global_labels_2d.ravel()
        global_labels_2d_flat[water_mask] = global_labels
        global_labels_2d = global_labels_2d_flat.reshape(H, W)
        
        return labels_2d, global_labels_2d, water_mask, H, W
    
    def find_scene_index(self, filename):
        for idx, scene_result in enumerate(self.scene_results):
            if scene_result.get('filename') == filename:
                return idx
        return None
    
    def update_display(self):
        if self.current_data is None or self.current_global_labels is None:
            return
        
        self.fig.clear()
        self.ax_left = self.fig.add_axes([0.02, 0.05, 0.58, 0.9])
        self.ax_right = self.fig.add_axes([0.65, 0.1, 0.32, 0.8])
        
        if self.show_flh:
            flh_masked = np.ma.masked_where(np.isnan(self.current_flh_2d), self.current_flh_2d)
            im_left = self.ax_left.imshow(flh_masked, cmap='jet',
                                         vmin=np.nanpercentile(self.current_flh, 2),
                                         vmax=np.nanpercentile(self.current_flh, 98))
            self.ax_left.set_title('FLH', fontsize=9, pad=3)
            
            cbar_ax = self.fig.add_axes([0.01, 0.15, 0.008, 0.7])
            cbar_left = self.fig.colorbar(im_left, cax=cbar_ax)
            cbar_left.set_label('FLH', fontsize=7, rotation=90, labelpad=8)
            cbar_left.ax.tick_params(labelsize=6)
        else:
            global_flh_masked = np.ma.masked_where(np.isnan(self.current_global_flh_map), self.current_global_flh_map)
            
            im_left = self.ax_left.imshow(global_flh_masked, cmap='RdYlBu_r',
                                         vmin=np.nanpercentile(self.current_global_flh_map[~np.isnan(self.current_global_flh_map)], 5),
                                         vmax=np.nanpercentile(self.current_global_flh_map[~np.isnan(self.current_global_flh_map)], 95))
            
            if len(self.highlighted_clusters) > 0:
                H, W = self.current_global_labels.shape
                highlight_mask = np.zeros((H, W), dtype=bool)
                for cluster_id in self.highlighted_clusters:
                    highlight_mask |= (self.current_global_labels == cluster_id)
                
                highlight_rgba = np.zeros((H, W, 4))
                highlight_rgba[highlight_mask] = [1.0, 1.0, 0.0, 1.0]
                self.ax_left.imshow(highlight_rgba)
            
            edges_h = sobel(self.current_global_labels, axis=0, mode='constant')
            edges_v = sobel(self.current_global_labels, axis=1, mode='constant')
            edges = np.hypot(edges_h, edges_v)
            edges = (edges > 0) & (self.current_global_labels >= 0)
            self.ax_left.contour(edges, levels=[0.5], colors='black', linewidths=0.5, alpha=0.8)
            
            self.ax_left.set_title('Global Clusters by Mean FLH', fontsize=9, pad=3)
            
            cbar_ax = self.fig.add_axes([0.01, 0.15, 0.008, 0.7])
            cbar_left = self.fig.colorbar(im_left, cax=cbar_ax)
            cbar_left.set_label('FLH', fontsize=7, rotation=90, labelpad=8)
            cbar_left.ax.tick_params(labelsize=6)
        
        self.ax_left.axis('off')
        
        self.ax_right.set_title('Rrs Spectra', fontsize=9, pad=3)
        self.ax_right.set_xlabel('Wavelength (nm)', fontsize=8)
        self.ax_right.set_ylabel('Rrs (sr^-1)', fontsize=8)
        self.ax_right.grid(True, alpha=0.3)
        self.ax_right.tick_params(labelsize=7)
        
        if len(self.spectra_list) > 0:
            for cluster_id, mean_rrs in self.spectra_list:
                self.ax_right.plot(self.wavelengths, mean_rrs, marker='o', 
                                  label=f'C{cluster_id}', linewidth=1.5, markersize=4)
            self.ax_right.legend(loc='best', fontsize=7, framealpha=0.9)
        else:
            self.ax_right.text(0.5, 0.5, 'Click on a cluster to view spectrum', 
                             ha='center', va='center', transform=self.ax_right.transAxes, fontsize=10)
        
        self.plot_nfis_points()
        self.canvas.draw()
    
    def on_mouse_press(self, event):
        if event.inaxes == self.ax_left and not self.show_flh:
            if event.xdata is not None and event.ydata is not None:
                x, y = int(event.xdata), int(event.ydata)
                H, W = self.current_global_labels.shape
                if 0 <= y < H and 0 <= x < W:
                    if self.current_global_labels[y, x] < 0:
                        self.show_flh = True
                        self.update_display()
    
    def on_mouse_release(self, event):
        if self.show_flh:
            self.show_flh = False
            self.update_display()
    
    def on_click(self, event):
        if event.inaxes == self.ax_left and event.xdata is not None and event.ydata is not None:
            x, y = int(event.xdata), int(event.ydata)
            H, W = self.current_global_labels.shape
            if 0 <= y < H and 0 <= x < W:
                cluster_id = self.current_global_labels[y, x]
                if cluster_id >= 0:
                    if self.delete_mode:
                        self.delete_cluster_from_archive(cluster_id)
                        self.delete_mode = False
                        self.delete_cluster_button.config(text="Delete Cluster")
                        return
                    
                    global_labels_flat = self.current_global_labels.ravel()[self.water_mask]
                    mask_cluster = (global_labels_flat == cluster_id)
                    if mask_cluster.sum() > 0:
                        mean_rrs = self.current_Rrs_all[mask_cluster].mean(axis=0)
                        
                        if not self.hold_var.get():
                            self.spectra_list = [(cluster_id, mean_rrs)]
                            self.highlighted_clusters = {cluster_id}
                        else:
                            existing_ids = [cid for cid, _ in self.spectra_list]
                            if cluster_id not in existing_ids:
                                self.spectra_list.append((cluster_id, mean_rrs))
                                self.highlighted_clusters.add(cluster_id)
                        
                        if len(self.highlighted_clusters) > 0:
                            self.archive_button.config(state=tk.NORMAL)
                        
                        self.update_display()
    
    def on_clear_spectra(self):
        self.spectra_list = []
        self.highlighted_clusters = set()
        self.archive_button.config(state=tk.DISABLED)
        self.update_display()
    
    def display_archived_clusters(self):
        if self.current_data is None or self.current_global_labels is None:
            messagebox.showwarning("Warning", "No clustering data available")
            return
        
        try:
            date_str = self.date_combo.get()
            archive_file = self.cluster_archive_dir / f'clusters_{date_str}.npz'
            
            if not archive_file.exists():
                messagebox.showinfo("Info", f"No archived clusters found for {date_str}")
                return
            
            existing_data = np.load(archive_file, allow_pickle=True)
            archived_clusters = existing_data['archived_clusters'].item()
            
            if len(archived_clusters) == 0:
                messagebox.showinfo("Info", f"No archived clusters found for {date_str}")
                return
            
            global_labels_flat = self.current_global_labels.ravel()[self.water_mask]
            unique_global_in_scene = set(np.unique(global_labels_flat[global_labels_flat >= 0]))
            
            archived_ids_in_scene = []
            for cluster_id in archived_clusters.keys():
                if cluster_id in unique_global_in_scene:
                    archived_ids_in_scene.append(cluster_id)
            
            if len(archived_ids_in_scene) == 0:
                messagebox.showinfo("Info", f"No archived clusters visible in current scene")
                return
            
            self.highlighted_clusters = set(archived_ids_in_scene)
            
            self.spectra_list = []
            for cluster_id in archived_ids_in_scene:
                mask_cluster = (global_labels_flat == cluster_id)
                if mask_cluster.sum() > 0:
                    mean_rrs = self.current_Rrs_all[mask_cluster].mean(axis=0)
                    self.spectra_list.append((cluster_id, mean_rrs))
            
            self.archive_button.config(state=tk.NORMAL)
            self.update_display()
            
            msg = f"Displaying {len(archived_ids_in_scene)} archived clusters\n"
            msg += f"Total archived for {date_str}: {len(archived_clusters)}"
            messagebox.showinfo("Archived Clusters", msg)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to display archived clusters:\n{str(e)}")
            import traceback
            traceback.print_exc()
    
    def archive_selected_clusters(self):
        if len(self.highlighted_clusters) == 0:
            messagebox.showwarning("Warning", "No clusters selected")
            return
        
        if self.current_data is None or self.current_global_labels is None:
            messagebox.showwarning("Warning", "No clustering data available")
            return
        
        try:
            date_str = self.date_combo.get()
            archive_file = self.cluster_archive_dir / f'clusters_{date_str}.npz'
            
            if archive_file.exists():
                existing_data = np.load(archive_file, allow_pickle=True)
                archived_clusters = existing_data['archived_clusters'].item()
            else:
                archived_clusters = {}
            
            global_labels_flat = self.current_global_labels.ravel()[self.water_mask]
            
            new_clusters_added = 0
            for cluster_id in self.highlighted_clusters:
                if cluster_id in archived_clusters:
                    print(f"Cluster {cluster_id} already archived, skipping...")
                    continue
                
                mask_cluster = (global_labels_flat == cluster_id)
                if mask_cluster.sum() == 0:
                    continue
                
                H, W = self.current_global_labels.shape
                
                water_indices = np.where(self.water_mask)[0]
                cluster_indices = water_indices[mask_cluster]
                
                rows = cluster_indices // W
                cols = cluster_indices % W
                
                lat_all = self.current_data['latitude'].ravel()
                lon_all = self.current_data['longitude'].ravel()
                
                lats = lat_all[cluster_indices]
                lons = lon_all[cluster_indices]
                
                rrs_data = {}
                for wl in self.wavelengths:
                    band_all = self.current_data[f'Rrs_{wl}'].ravel()
                    rrs_data[f'Rrs_{wl}'] = band_all[cluster_indices]
                
                rhoc_data = {}
                for wl in self.wavelengths:
                    band_all = self.current_data[f'RhoC_{wl}'].ravel()
                    rhoc_data[f'RhoC_{wl}'] = band_all[cluster_indices]
                
                band_660 = self.current_data['Rrs_660'].ravel()[cluster_indices]
                band_709 = self.current_data['Rrs_709'].ravel()[cluster_indices]
                band_745 = self.current_data['Rrs_745'].ravel()[cluster_indices]
                flh = compute_flh(band_660, band_709, band_745)
                
                cluster_info = {
                    'global_cluster_id': int(cluster_id),
                    'n_pixels': int(mask_cluster.sum()),
                    'latitude': lats,
                    'longitude': lons,
                    'rows': rows,
                    'cols': cols,
                    'source_file': self.current_file.name,
                    'date': date_str,
                    'timestamp': datetime.now().strftime('%Y%m%d_%H%M%S')
                }
                
                cluster_info.update(rrs_data)
                cluster_info.update(rhoc_data)
                cluster_info['FLH'] = flh
                
                archived_clusters[cluster_id] = cluster_info
                new_clusters_added += 1
            
            np.savez_compressed(archive_file, archived_clusters=archived_clusters)
            
            self.generate_archive_summary(date_str, archived_clusters)
            
            total_clusters = len(archived_clusters)
            msg = f"Archive updated!\n"
            msg += f"New clusters added: {new_clusters_added}\n"
            msg += f"Total clusters in {date_str}: {total_clusters}\n"
            msg += f"File: {archive_file.name}"
            
            messagebox.showinfo("Success", msg)
            print(f"Archived {new_clusters_added} clusters to {archive_file}")
            
            self.display_archived_button.config(state=tk.NORMAL)
            self.delete_cluster_button.config(state=tk.NORMAL)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to archive clusters:\n{str(e)}")
            import traceback
            traceback.print_exc()
    
    def toggle_delete_mode(self):
        self.delete_mode = not self.delete_mode
        if self.delete_mode:
            self.delete_cluster_button.config(text="Delete Cluster (Active)")
            self.info_label.config(text="Click on a cluster to delete from archive")
        else:
            self.delete_cluster_button.config(text="Delete Cluster")
            self.info_label.config(text="Delete mode cancelled")
    
    def delete_cluster_from_archive(self, cluster_id):
        try:
            date_str = self.date_combo.get()
            archive_file = self.cluster_archive_dir / f'clusters_{date_str}.npz'
            
            if not archive_file.exists():
                messagebox.showwarning("Warning", f"No archive file found for {date_str}")
                return
            
            existing_data = np.load(archive_file, allow_pickle=True)
            archived_clusters = existing_data['archived_clusters'].item()
            
            if cluster_id not in archived_clusters:
                messagebox.showinfo("Info", f"Cluster {cluster_id} is not in archive")
                return
            
            cluster_info = archived_clusters[cluster_id]
            result = messagebox.askyesno("Confirm Delete", 
                                        f"Delete cluster {cluster_id} from archive?\n"
                                        f"Source: {cluster_info.get('source_file', 'unknown')}\n"
                                        f"Pixels: {cluster_info.get('n_pixels', 0)}")
            
            if result:
                del archived_clusters[cluster_id]
                
                np.savez_compressed(archive_file, archived_clusters=archived_clusters)
                
                if cluster_id in self.highlighted_clusters:
                    self.highlighted_clusters.remove(cluster_id)
                    self.spectra_list = [(cid, rrs) for cid, rrs in self.spectra_list if cid != cluster_id]
                    self.update_display()
                
                msg = f"Cluster {cluster_id} deleted from archive\n"
                msg += f"Remaining clusters: {len(archived_clusters)}"
                messagebox.showinfo("Success", msg)
                print(f"Deleted cluster {cluster_id} from {archive_file}")
                
                if len(archived_clusters) == 0:
                    self.delete_cluster_button.config(state=tk.DISABLED)
                else:
                    date_str = self.date_combo.get()
                    self.generate_archive_summary(date_str, archived_clusters)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to delete cluster:\n{str(e)}")
            import traceback
            traceback.print_exc()
    
    def generate_archive_summary(self, date_str, archived_clusters):
        try:
            summary_file = self.cluster_archive_dir / f'clusters_{date_str}.txt'
            
            by_source = {}
            total_pixels = 0
            
            for cluster_id, cluster_info in archived_clusters.items():
                source_file = cluster_info.get('source_file', 'unknown')
                n_pixels = cluster_info.get('n_pixels', 0)
                
                if source_file not in by_source:
                    by_source[source_file] = {
                        'cluster_ids': [],
                        'total_pixels': 0
                    }
                
                by_source[source_file]['cluster_ids'].append(cluster_id)
                by_source[source_file]['total_pixels'] += n_pixels
                total_pixels += n_pixels
            
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write(f"Archive Summary for Date: {date_str}\n")
                f.write(f"{'='*70}\n\n")
                
                f.write(f"Total Clusters: {len(archived_clusters)}\n")
                f.write(f"Total Pixels: {total_pixels:,}\n")
                f.write(f"Number of Source Files: {len(by_source)}\n\n")
                
                f.write(f"{'='*70}\n")
                f.write(f"Breakdown by Source File:\n")
                f.write(f"{'='*70}\n\n")
                
                for source_file in sorted(by_source.keys()):
                    info = by_source[source_file]
                    f.write(f"Source: {source_file}\n")
                    f.write(f"  Number of Clusters: {len(info['cluster_ids'])}\n")
                    f.write(f"  Total Pixels: {info['total_pixels']:,}\n")
                    f.write(f"  Cluster IDs: {sorted(info['cluster_ids'])}\n")
                    f.write(f"\n")
                
                f.write(f"{'='*70}\n")
                f.write(f"Detailed Cluster Information:\n")
                f.write(f"{'='*70}\n\n")
                
                for cluster_id in sorted(archived_clusters.keys()):
                    cluster_info = archived_clusters[cluster_id]
                    f.write(f"Cluster ID: {cluster_id}\n")
                    f.write(f"  Source File: {cluster_info.get('source_file', 'unknown')}\n")
                    f.write(f"  Number of Pixels: {cluster_info.get('n_pixels', 0):,}\n")
                    f.write(f"  Timestamp: {cluster_info.get('timestamp', 'unknown')}\n")
                    
                    if 'latitude' in cluster_info and len(cluster_info['latitude']) > 0:
                        lats = cluster_info['latitude']
                        lons = cluster_info['longitude']
                        f.write(f"  Latitude Range: [{lats.min():.4f}, {lats.max():.4f}]\n")
                        f.write(f"  Longitude Range: [{lons.min():.4f}, {lons.max():.4f}]\n")
                    
                    f.write(f"\n")
            
            print(f"Generated summary: {summary_file}")
            
        except Exception as e:
            print(f"Warning: Failed to generate summary: {e}")
            import traceback
            traceback.print_exc()


if __name__ == '__main__':
    root = tk.Tk()
    app = ClusteringGUI(root)
    root.mainloop()

