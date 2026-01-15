"""
P01_GUI_region_selector.py

GOCI-II Region Selector GUI
- Display GOCI-II RhoC RGB
- Draw rectangular regions
- Save regions to JSON
- Crop regions to NPZ files
- Load and display NFIS points
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from pathlib import Path
from datetime import datetime
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import json
import pandas as pd

# Add library path
sys.path.insert(0, str(Path(__file__).parent / 'library'))

from goci2_reader import read_goci2


class RegionSelector:
    def __init__(self, root):
        self.root = root
        self.root.title("GOCI-II Region Selector")
        self.root.geometry("1600x900")
        
        self.data = None
        self.current_rect = None
        self.press = None
        self.current_selection = None
        self.saved_regions = []
        self.displayed_rect = None
        self.crop_offset = None
        self.nfis_scatter = None
        self.nfis_data = None
        
        # Paths
        self.base_dir = Path(__file__).parent
        self.data_dir = self.base_dir / 'data' / 'goci2'
        self.nfis_dir = self.base_dir / 'data' / 'nfis'
        self.results_dir = self.base_dir / 'results' / 'GUI_archive'
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.regions_file = self.results_dir / 'regions.json'
        
        # Find available GOCI-II files
        self.available_files = sorted([f.name for f in self.data_dir.glob('*.nc')])
        self.default_file = self.available_files[0] if self.available_files else None
        
        self.create_widgets()
        if self.default_file:
            self.load_data()
        self.load_saved_regions()
        
    def create_widgets(self):
        main_container = ttk.Frame(self.root)
        main_container.pack(fill=tk.BOTH, expand=True)
        
        # Left panel - canvas
        canvas_frame = ttk.Frame(main_container)
        canvas_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # File selector
        control_top = ttk.Frame(canvas_frame)
        control_top.pack(side=tk.TOP, fill=tk.X, pady=(0, 5))
        
        ttk.Label(control_top, text="Background Image:").pack(side=tk.LEFT, padx=(0, 5))
        self.file_selector = ttk.Combobox(control_top, values=self.available_files, 
                                          state='readonly', width=50)
        self.file_selector.pack(side=tk.LEFT, fill=tk.X, expand=True)
        if self.default_file:
            self.file_selector.set(self.default_file)
        self.file_selector.bind('<<ComboboxSelected>>', self.on_file_change)
        
        # Matplotlib canvas
        self.fig = Figure(figsize=(10, 8))
        self.ax = self.fig.add_subplot(111)
        self.ax.axis('off')
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        
        # Right panel - controls
        right_panel = ttk.Frame(main_container, width=350)
        right_panel.pack(side=tk.RIGHT, fill=tk.Y, padx=10, pady=10)
        right_panel.pack_propagate(False)
        
        # Info display
        info_frame = ttk.LabelFrame(right_panel, text="Current Rectangle Info", padding=10)
        info_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.info_text = tk.Text(info_frame, height=12, width=40, wrap=tk.WORD, state=tk.DISABLED)
        self.info_text.pack(fill=tk.BOTH, expand=True)
        
        # Save/Crop buttons
        save_frame = ttk.LabelFrame(right_panel, text="Region Actions", padding=10)
        save_frame.pack(fill=tk.X, pady=(0, 10))
        
        button_row = ttk.Frame(save_frame)
        button_row.pack(fill=tk.X, pady=(0, 5))
        
        self.save_button = ttk.Button(button_row, text="Save", command=self.show_save_dialog, state=tk.DISABLED)
        self.save_button.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 2))
        
        self.crop_button = ttk.Button(button_row, text="Crop", command=self.crop_current_region, state=tk.DISABLED)
        self.crop_button.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(2, 0))
        
        # Save dialog (hidden by default)
        self.save_dialog_frame = ttk.Frame(save_frame)
        
        ttk.Label(self.save_dialog_frame, text="Region Name:").pack(anchor=tk.W, pady=(10, 2))
        self.region_name_entry = ttk.Entry(self.save_dialog_frame, width=30)
        self.region_name_entry.pack(fill=tk.X, pady=(0, 5))
        
        button_frame = ttk.Frame(self.save_dialog_frame)
        button_frame.pack(fill=tk.X)
        ttk.Button(button_frame, text="Submit", command=self.submit_region).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 2))
        ttk.Button(button_frame, text="Cancel", command=self.cancel_save).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(2, 0))
        
        # Regions list
        regions_frame = ttk.LabelFrame(right_panel, text="Saved Regions", padding=10)
        regions_frame.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = ttk.Scrollbar(regions_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.regions_listbox = tk.Listbox(regions_frame, yscrollcommand=scrollbar.set)
        self.regions_listbox.pack(fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.regions_listbox.yview)
        
        self.regions_listbox.bind('<<ListboxSelect>>', self.on_region_select)
        
        delete_button = ttk.Button(regions_frame, text="Delete Selected", command=self.delete_selected_region)
        delete_button.pack(fill=tk.X, pady=(5, 0))
        
        # NFIS controls
        nfis_frame = ttk.LabelFrame(right_panel, text="NFIS Data", padding=10)
        nfis_frame.pack(fill=tk.X, pady=(10, 0))
        
        ttk.Button(nfis_frame, text="Add NFIS", command=self.load_nfis_data).pack(fill=tk.X, pady=(0, 5))
        ttk.Button(nfis_frame, text="Clear NFIS", command=self.clear_nfis_data).pack(fill=tk.X)
        
    def load_data(self, filename=None):
        if filename is None:
            filename = self.default_file
        
        target_file = self.data_dir / filename
        
        if not target_file.exists():
            messagebox.showerror("Error", f"File not found:\n{target_file}")
            return
        
        try:
            self.data = read_goci2(str(target_file))
            self.current_file = filename
            self.display_rrs()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")
    
    def on_file_change(self, event):
        selected_file = self.file_selector.get()
        if selected_file:
            if self.current_rect is not None:
                self.current_rect.remove()
                self.current_rect = None
            if self.displayed_rect is not None:
                self.displayed_rect.remove()
                self.displayed_rect = None
            self.current_selection = None
            self.save_button.config(state=tk.DISABLED)
            self.crop_button.config(state=tk.DISABLED)
            self.load_data(selected_file)
    
    def display_rrs(self):
        if self.data is None:
            return
        
        # Get RGB bands (RhoC)
        r_band = self.data.get_band('rhoc', 660)
        g_band = self.data.get_band('rhoc', 555)
        b_band = self.data.get_band('rhoc', 443)
        
        if r_band is None or g_band is None or b_band is None:
            messagebox.showerror("Error", "Failed to get RhoC bands")
            return
        
        # Crop to top 2/3 and right 2/3
        height, width = self.data.dimensions
        y_end = int(height * 2 / 3)
        x_start = int(width * 1 / 3)
        
        r_band = r_band[:y_end, x_start:]
        g_band = g_band[:y_end, x_start:]
        b_band = b_band[:y_end, x_start:]
        
        self.crop_offset = (0, y_end, x_start, width)
        
        # Normalize
        r_norm = self.normalize_band(r_band, 0.0, 0.2)
        g_norm = self.normalize_band(g_band, 0.0, 0.2)
        b_norm = self.normalize_band(b_band, 0.0, 0.2)
        
        rgb = np.dstack([r_norm, g_norm, b_norm])
        
        # Display
        self.ax.clear()
        self.ax.imshow(rgb)
        self.ax.axis('off')
        
        date_str = self.current_file[17:25] if hasattr(self, 'current_file') else 'GOCI-II'
        self.ax.set_title(f'GOCI-II RhoC RGB ({date_str})', fontsize=12, fontweight='bold')
        
        # Redisplay NFIS points if loaded
        if self.nfis_data is not None:
            self.nfis_scatter = None
            self.display_nfis_points()
        
        self.canvas.draw()
    
    def normalize_band(self, band, vmin, vmax):
        valid_mask = np.isfinite(band)
        normalized = np.clip((band - vmin) / (vmax - vmin), 0, 1)
        normalized[~valid_mask] = 0
        return normalized
    
    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        
        # Clear existing rectangles
        if self.current_rect is not None:
            self.current_rect.remove()
            self.current_rect = None
        
        if self.displayed_rect is not None:
            self.displayed_rect.remove()
            self.displayed_rect = None
        
        self.press = (event.xdata, event.ydata)
        self.canvas.draw()
    
    def on_motion(self, event):
        if self.press is None or event.inaxes != self.ax:
            return
        
        if self.current_rect is not None:
            self.current_rect.remove()
        
        x0, y0 = self.press
        x1, y1 = event.xdata, event.ydata
        
        width = x1 - x0
        height = y1 - y0
        
        self.current_rect = Rectangle((x0, y0), width, height,
                                     linewidth=2, edgecolor='yellow',
                                     facecolor='none', linestyle='-')
        self.ax.add_patch(self.current_rect)
        self.canvas.draw()
    
    def on_release(self, event):
        if self.press is None or event.inaxes != self.ax:
            return
        
        x0, y0 = self.press
        x1, y1 = event.xdata, event.ydata
        self.press = None
        
        # Ignore too small rectangles
        if abs(x1 - x0) < 5 or abs(y1 - y0) < 5:
            if self.current_rect is not None:
                self.current_rect.remove()
                self.current_rect = None
                self.canvas.draw()
            return
        
        x_min, x_max = min(x0, x1), max(x0, x1)
        y_min, y_max = min(y0, y1), max(y0, y1)
        
        # Convert to pixel coordinates
        col_start = int(np.round(x_min))
        col_end = int(np.round(x_max))
        row_start = int(np.round(y_min))
        row_end = int(np.round(y_max))
        
        # Convert to original image coordinates
        if self.crop_offset is not None:
            y_offset_start, y_offset_end, x_offset_start, x_offset_end = self.crop_offset
            col_start += x_offset_start
            col_end += x_offset_start
            row_start += y_offset_start
            row_end += y_offset_start
        
        # Clamp to valid range
        height, width = self.data.dimensions
        col_start = max(0, min(col_start, width - 1))
        col_end = max(0, min(col_end, width))
        row_start = max(0, min(row_start, height - 1))
        row_end = max(0, min(row_end, height))
        
        self.current_selection = {
            'row_start': row_start,
            'row_end': row_end,
            'col_start': col_start,
            'col_end': col_end
        }
        
        self.update_info_display()
        self.save_button.config(state=tk.NORMAL)
        self.crop_button.config(state=tk.NORMAL)
    
    def update_info_display(self):
        if self.current_selection is None:
            return
        
        row_start = self.current_selection['row_start']
        row_end = self.current_selection['row_end']
        col_start = self.current_selection['col_start']
        col_end = self.current_selection['col_end']
        
        lat_top_left = self.data.latitude[row_start, col_start]
        lon_top_left = self.data.longitude[row_start, col_start]
        
        width_pixels = col_end - col_start
        height_pixels = row_end - row_start
        total_pixels = width_pixels * height_pixels
        
        info = f"Top-left Latitude: {lat_top_left:.6f}\n"
        info += f"Top-left Longitude: {lon_top_left:.6f}\n\n"
        info += f"Top-left Row: {row_start}\n"
        info += f"Top-left Column: {col_start}\n\n"
        info += f"Width (pixels): {width_pixels}\n"
        info += f"Height (pixels): {height_pixels}\n\n"
        info += f"Total Pixels: {total_pixels}"
        
        self.info_text.config(state=tk.NORMAL)
        self.info_text.delete(1.0, tk.END)
        self.info_text.insert(1.0, info)
        self.info_text.config(state=tk.DISABLED)
    
    def show_save_dialog(self):
        if self.current_selection is None:
            return
        
        self.save_dialog_frame.pack(fill=tk.X, pady=(10, 0))
        self.region_name_entry.focus()
    
    def cancel_save(self):
        self.save_dialog_frame.pack_forget()
        self.region_name_entry.delete(0, tk.END)
    
    def submit_region(self):
        region_name = self.region_name_entry.get().strip()
        
        if not region_name:
            messagebox.showwarning("Warning", "Please enter a region name")
            return
        
        if any(r['name'] == region_name for r in self.saved_regions):
            messagebox.showwarning("Warning", f"Region name '{region_name}' already exists")
            return
        
        if self.current_selection is None:
            return
        
        row_start = self.current_selection['row_start']
        row_end = self.current_selection['row_end']
        col_start = self.current_selection['col_start']
        col_end = self.current_selection['col_end']
        
        lat_top_left = float(self.data.latitude[row_start, col_start])
        lon_top_left = float(self.data.longitude[row_start, col_start])
        
        region_data = {
            'name': region_name,
            'row_start': row_start,
            'row_end': row_end,
            'col_start': col_start,
            'col_end': col_end,
            'lat_top_left': lat_top_left,
            'lon_top_left': lon_top_left,
            'width_pixels': col_end - col_start,
            'height_pixels': row_end - row_start,
            'total_pixels': (col_end - col_start) * (row_end - row_start),
            'timestamp': datetime.now().strftime('%Y%m%d_%H%M%S')
        }
        
        self.saved_regions.append(region_data)
        self.save_regions_to_file()
        self.update_regions_list()
        
        self.cancel_save()
        messagebox.showinfo("Success", f"Region '{region_name}' saved successfully")
    
    def save_regions_to_file(self):
        with open(self.regions_file, 'w') as f:
            json.dump(self.saved_regions, f, indent=2)
    
    def load_saved_regions(self):
        if self.regions_file.exists():
            try:
                with open(self.regions_file, 'r') as f:
                    self.saved_regions = json.load(f)
                self.update_regions_list()
            except Exception as e:
                print(f"Failed to load saved regions: {e}")
    
    def update_regions_list(self):
        self.regions_listbox.delete(0, tk.END)
        for region in self.saved_regions:
            self.regions_listbox.insert(tk.END, region['name'])
    
    def on_region_select(self, event):
        selection = self.regions_listbox.curselection()
        if not selection:
            return
        
        index = selection[0]
        region = self.saved_regions[index]
        
        # Clear current rectangle
        if self.current_rect is not None:
            self.current_rect.remove()
            self.current_rect = None
        
        # Remove old displayed rectangle
        if self.displayed_rect is not None:
            self.displayed_rect.remove()
        
        row_start = region['row_start']
        row_end = region['row_end']
        col_start = region['col_start']
        col_end = region['col_end']
        
        # Convert to display coordinates
        display_row_start = row_start
        display_row_end = row_end
        display_col_start = col_start
        display_col_end = col_end
        
        if self.crop_offset is not None:
            y_offset_start, y_offset_end, x_offset_start, x_offset_end = self.crop_offset
            display_row_start -= y_offset_start
            display_row_end -= y_offset_start
            display_col_start -= x_offset_start
            display_col_end -= x_offset_start
        
        width = display_col_end - display_col_start
        height = display_row_end - display_row_start
        
        self.displayed_rect = Rectangle((display_col_start, display_row_start), width, height,
                                       linewidth=2, edgecolor='cyan',
                                       facecolor='none', linestyle='-')
        self.ax.add_patch(self.displayed_rect)
        self.canvas.draw()
        
        self.current_selection = {
            'row_start': row_start,
            'row_end': row_end,
            'col_start': col_start,
            'col_end': col_end
        }
        self.update_info_display()
        self.save_button.config(state=tk.NORMAL)
        self.crop_button.config(state=tk.NORMAL)
    
    def delete_selected_region(self):
        selection = self.regions_listbox.curselection()
        if not selection:
            messagebox.showwarning("Warning", "Please select a region to delete")
            return
        
        index = selection[0]
        region_name = self.saved_regions[index]['name']
        
        if messagebox.askyesno("Confirm Delete", f"Delete region '{region_name}'?"):
            del self.saved_regions[index]
            self.save_regions_to_file()
            self.update_regions_list()
            
            if self.displayed_rect is not None:
                self.displayed_rect.remove()
                self.displayed_rect = None
                self.canvas.draw()
    
    def crop_current_region(self):
        if self.current_selection is None:
            messagebox.showwarning("Warning", "No region selected")
            return
        
        if not hasattr(self, 'current_file') or self.data is None:
            messagebox.showwarning("Warning", "No data loaded")
            return
        
        row_start = self.current_selection['row_start']
        row_end = self.current_selection['row_end']
        col_start = self.current_selection['col_start']
        col_end = self.current_selection['col_end']
        
        # Find region name
        region_name = None
        for region in self.saved_regions:
            if (region['row_start'] == row_start and region['row_end'] == row_end and
                region['col_start'] == col_start and region['col_end'] == col_end):
                region_name = region['name']
                break
        
        if region_name is None:
            result = messagebox.askquestion("Region Name", 
                                           "This region is not saved yet.\nDo you want to crop without a region name?",
                                           icon='question')
            if result == 'no':
                return
            region_name = "unnamed"
        
        # Crop all bands
        save_dict = {}
        
        wavelengths = [380, 412, 443, 490, 510, 555, 620, 660, 680, 709, 745, 865]
        
        for wl in wavelengths:
            rhoc_band = self.data.get_band('rhoc', wl)
            if rhoc_band is not None:
                cropped = rhoc_band[row_start:row_end, col_start:col_end]
                save_dict[f'RhoC_{wl}'] = cropped
            
            rrs_band = self.data.get_band('rrs', wl)
            if rrs_band is not None:
                cropped = rrs_band[row_start:row_end, col_start:col_end]
                save_dict[f'Rrs_{wl}'] = cropped
        
        # Add geolocation and flag
        save_dict['latitude'] = self.data.latitude[row_start:row_end, col_start:col_end]
        save_dict['longitude'] = self.data.longitude[row_start:row_end, col_start:col_end]
        save_dict['flag'] = self.data.flag[row_start:row_end, col_start:col_end]
        
        # Calculate FLH
        flh_rhoc = self.calculate_flh('rhoc')
        if flh_rhoc is not None:
            save_dict['FLH_RhoC'] = flh_rhoc[row_start:row_end, col_start:col_end]
        
        flh_rrs = self.calculate_flh('rrs')
        if flh_rrs is not None:
            save_dict['FLH_Rrs'] = flh_rrs[row_start:row_end, col_start:col_end]
        
        # Save to file
        crop_dir = self.results_dir / 'cropped_data'
        crop_dir.mkdir(parents=True, exist_ok=True)
        
        base_name = Path(self.current_file).stem
        output_filename = f"{base_name}_cropped_{row_start}-{row_end}_{col_start}-{col_end}_{region_name}.npz"
        output_path = crop_dir / output_filename
        
        np.savez_compressed(output_path, **save_dict)
        
        info_msg = f"Cropped data saved to:\n{output_path.name}\n\n"
        info_msg += f"Region: Y[{row_start}:{row_end}], X[{col_start}:{col_end}]\n"
        info_msg += f"Shape: {row_end-row_start} x {col_end-col_start}\n"
        info_msg += f"RhoC bands: {len([k for k in save_dict.keys() if k.startswith('RhoC')])}\n"
        info_msg += f"Rrs bands: {len([k for k in save_dict.keys() if k.startswith('Rrs')])}\n"
        
        messagebox.showinfo("Success", info_msg)
    
    def calculate_flh(self, data_type):
        """Calculate FLH with baseline 660-745nm"""
        band_660 = self.data.get_band(data_type, 660)
        band_709 = self.data.get_band(data_type, 709)
        band_745 = self.data.get_band(data_type, 745)
        
        if band_660 is None or band_709 is None or band_745 is None:
            return None
        
        lambda_660 = 660.0
        lambda_709 = 709.0
        lambda_745 = 745.0
        
        baseline = band_660 + (band_745 - band_660) * (lambda_709 - lambda_660) / (lambda_745 - lambda_660)
        flh = band_709 - baseline
        
        return flh
    
    def load_nfis_data(self):
        filepath = filedialog.askopenfilename(
            title="Select NFIS file",
            initialdir=str(self.nfis_dir),
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        
        if not filepath:
            return
        
        try:
            if filepath.endswith('.csv'):
                df = pd.read_csv(filepath)
            else:
                df = pd.read_excel(filepath)
            
            if 'latitude' not in df.columns or 'longitude' not in df.columns:
                messagebox.showerror("Error", "File must contain 'latitude' and 'longitude' columns")
                return
            
            self.nfis_data = df[['latitude', 'longitude']].dropna()
            self.display_nfis_points()
            messagebox.showinfo("Success", f"Loaded {len(self.nfis_data)} NFIS points")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load NFIS file:\n{str(e)}")
    
    def clear_nfis_data(self):
        if self.nfis_scatter is not None:
            self.nfis_scatter.remove()
            self.nfis_scatter = None
            self.nfis_data = None
            self.canvas.draw()
    
    def display_nfis_points(self):
        if self.nfis_data is None or self.data is None:
            return
        
        if self.nfis_scatter is not None:
            self.nfis_scatter.remove()
            self.nfis_scatter = None
        
        lat_points = self.nfis_data['latitude'].values
        lon_points = self.nfis_data['longitude'].values
        
        lat_grid = self.data.latitude
        lon_grid = self.data.longitude
        
        # Find nearest pixels
        pixel_rows = []
        pixel_cols = []
        
        for lat, lon in zip(lat_points, lon_points):
            distances = np.sqrt((lat_grid - lat)**2 + (lon_grid - lon)**2)
            min_idx = np.unravel_index(np.argmin(distances), distances.shape)
            pixel_rows.append(min_idx[0])
            pixel_cols.append(min_idx[1])
        
        # Convert to display coordinates
        display_rows = []
        display_cols = []
        
        if self.crop_offset is not None:
            y_offset_start, y_offset_end, x_offset_start, x_offset_end = self.crop_offset
            
            for row, col in zip(pixel_rows, pixel_cols):
                if y_offset_start <= row < y_offset_end and x_offset_start <= col < x_offset_end:
                    display_rows.append(row - y_offset_start)
                    display_cols.append(col - x_offset_start)
        else:
            display_rows = pixel_rows
            display_cols = pixel_cols
        
        # Plot
        if display_rows:
            self.nfis_scatter = self.ax.scatter(display_cols, display_rows, 
                                               c='red', s=50, marker='o', 
                                               alpha=0.7, edgecolors='darkred', linewidths=1)
            self.canvas.draw()


if __name__ == '__main__':
    root = tk.Tk()
    app = RegionSelector(root)
    root.mainloop()

