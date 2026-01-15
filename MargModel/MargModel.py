# -*- coding: utf-8 -*-
"""
MargModel.py - Standalone Forward Model for CCRR-RT

This module provides a standalone implementation of the CCRR-RT forward model
used in the K03 optimization. It can be distributed independently to allow
third parties to compute Rrs spectra given water quality parameters.

Usage Example:
    import numpy as np
    import MargModel
    
    # Define wavelengths
    wavelengths = np.arange(400., 750.)
    
    # Set water quality parameters
    chla = 2.0      # Chlorophyll-a concentration (mg/m3)
    mineral = 0.5   # Mineral concentration (g/m3)
    aCDOM = 0.02    # CDOM absorption at 443nm (1/m)
    
    # Run forward model
    Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)
    
    # Rrs now contains the modeled remote sensing reflectance

Author: wonk1
Date: 2024
"""

import numpy as np
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import pickle
import os


# Global cache variables for performance optimization
_AB_BOX_CACHE = None
_AB_INTERP_CACHE = {}
_WAVELENGTH_EXP_CACHE = {}
_MOREL_INTERPOLATORS = None
_WATER_AB_CACHE = None


def _get_lamb_key(lamb):
    """
    Fast wavelength array key generation for caching.
    Uses shape + first/last values instead of full tuple conversion.
    """
    if isinstance(lamb, np.ndarray):
        return (lamb.shape[0], float(lamb[0]), float(lamb[-1]))
    else:
        return lamb


def interpolate_data(wl, Lsky, wl_ref):
    """
    Interpolate data using cubic spline.
    
    Parameters:
    -----------
    wl : array
        Original wavelengths
    Lsky : array
        Original data values
    wl_ref : array
        Target wavelengths for interpolation
        
    Returns:
    --------
    Lsky_new : array
        Interpolated data at wl_ref wavelengths
    """
    ind_good = ~(np.isnan(wl) | np.isnan(Lsky))
    tck = interpolate.splrep(wl[ind_good], Lsky[ind_good], s=0)
    Lsky_new = interpolate.splev(wl_ref, tck, der=0)
    return Lsky_new


def get_water_ab(lamb):
    """
    Get water absorption and backscattering coefficients.
    Cached to avoid repeated file I/O.
    
    Parameters:
    -----------
    lamb : array
        Wavelengths (nm)
        
    Returns:
    --------
    aw : array
        Water absorption coefficient (1/m)
    bw : array
        Water scattering coefficient (1/m)
    """
    global _WATER_AB_CACHE
    
    if _WATER_AB_CACHE is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        fname_ab = os.path.join(current_dir, "independent", "Water_a_b.txt")
        
        wl = []
        aa = []
        bb = []
        
        with open(fname_ab, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                words = line.split()
                if len(words) < 3:
                    continue
                wl.append(float(words[0]))
                aa.append(float(words[1]))
                bb.append(float(words[2]))
        
        wl_ab = np.array(wl)
        aw_old = np.array(aa)
        bw_old = np.array(bb)
        _WATER_AB_CACHE = (wl_ab, aw_old, bw_old)
    
    wl_ab, aw_old, bw_old = _WATER_AB_CACHE
    aw = interpolate_data(wl_ab, aw_old, lamb)
    bw = interpolate_data(wl_ab, bw_old, lamb)
    
    return aw, bw


# Morel 2002 bidirectional reflectance distribution function (BRDF) tables
wl_fq = np.array([412.5, 442.5, 490., 510., 560., 620., 660])
chla_fq = np.array([0.03, 0.1, 0.3, 1., 3., 10.])

fq_morel = np.array([
    [0.089538, 0.094417, 0.106584, 0.112755, 0.118970, 0.118219, 0.118756], 
    [0.095834, 0.097003, 0.101338, 0.104170, 0.107362, 0.105764, 0.104753],
    [0.097576, 0.096815, 0.097281, 0.097115, 0.097048, 0.094493, 0.093036], 
    [0.095803, 0.094307, 0.094060, 0.092853, 0.091167, 0.086045, 0.084662],
    [0.092191, 0.090745, 0.091792, 0.091698, 0.092447, 0.082218, 0.080325],
    [0.087741, 0.086690, 0.089566, 0.091422, 0.098468, 0.082820, 0.078632]])

sfq_morel = np.array([
    [0.001395, 0.002067, 0.003802, 0.005239, 0.005314, 0.005261, 0.005363],
    [0.000254, 0.000101, 0.000308, 0.000721, 0.000318, 0.000162, 0.000289],
    [0.001678, 0.002867, 0.005181, 0.004966, 0.005726, 0.004667, 0.005080],
    [0.004727, 0.006759, 0.010777, 0.011601, 0.013461, 0.009248, 0.008975],
    [0.008255, 0.010688, 0.015709, 0.017363, 0.020826, 0.015052, 0.013596],
    [0.012713, 0.015235, 0.021056, 0.023357, 0.028172, 0.023455, 0.020462]])

f0_morel = np.array([
    [0.297892, 0.311742, 0.347280, 0.359728, 0.375008, 0.370053, 0.372716],
    [0.324018, 0.328848, 0.345755, 0.350503, 0.358735, 0.350497, 0.349206],
    [0.340239, 0.341657, 0.350980, 0.349334, 0.349570, 0.334437, 0.330755],
    [0.351673, 0.352505, 0.362207, 0.359773, 0.357388, 0.327275, 0.320913],
    [0.359587, 0.360429, 0.374357, 0.376513, 0.383433, 0.335178, 0.323731],
    [0.370570, 0.370782, 0.389128, 0.397841, 0.424316, 0.362306, 0.34020]])

sf_morel = np.array([
    [0.065801, 0.076526, 0.095435, 0.103901, 0.121165, 0.134426, 0.143912],
    [0.095786, 0.111534, 0.138252, 0.147280, 0.169210, 0.183429, 0.195362],
    [0.131988, 0.154209, 0.191203, 0.201694, 0.225606, 0.227381, 0.236441],
    [0.183170, 0.213187, 0.261757, 0.276009, 0.305781, 0.284162, 0.285455],
    [0.239626, 0.273898, 0.330935, 0.349059, 0.386471, 0.352961, 0.343078],
    [0.316124, 0.351720, 0.415626, 0.438584, 0.482277, 0.457638, 0.433943]])

q0_morel = np.array([
    [3.318220, 3.291250, 3.245640, 3.176500, 3.138020, 3.116680, 3.126530],
    [3.375400, 3.385700, 3.408430, 3.359680, 3.336890, 3.309970, 3.332380],
    [3.484950, 3.529680, 3.613410, 3.601660, 3.606950, 3.542300, 3.560200],
    [3.675060, 3.746830, 3.868930, 3.894090, 3.942030, 3.815240, 3.801910],
    [3.913530, 3.991380, 4.110030, 4.141140, 4.188650, 4.104160, 4.053610],
    [4.252700, 4.313260, 4.393950, 4.405460, 4.368130, 4.427190, 4.371050]])

sq_morel = np.array([
    [0.863223, 0.976278, 1.129290, 1.203100, 1.296770, 1.413010, 1.479420],
    [1.055510, 1.190090, 1.382130, 1.482910, 1.625960, 1.775070, 1.866290],
    [1.302830, 1.469930, 1.700680, 1.827460, 2.036100, 2.175700, 2.269360],
    [1.671180, 1.877700, 2.115910, 2.239960, 2.479820, 2.711180, 2.789710],
    [2.083950, 2.303690, 2.506640, 2.580090, 2.707700, 3.141310, 3.229040],
    [2.625950, 2.843750, 2.978790, 2.986960, 2.898750, 3.527960, 3.717540]])


def get_morel_numbers(wl, chla, method='linear'):
    """
    Get Morel BRDF numbers with cached interpolators for efficiency.
    
    Parameters:
    -----------
    wl : float
        Wavelength (nm)
    chla : float
        Chlorophyll-a concentration (mg/m3)
    method : str
        Interpolation method ('linear' or 'cubic')
        
    Returns:
    --------
    res : list
        [f0, sf, q0, sq, fq, sfq] Morel parameters
    """
    global _MOREL_INTERPOLATORS
    
    if _MOREL_INTERPOLATORS is None:
        _MOREL_INTERPOLATORS = {
            'f0': RegularGridInterpolator((chla_fq, wl_fq), f0_morel, method=method),
            'sf': RegularGridInterpolator((chla_fq, wl_fq), sf_morel, method=method),
            'q0': RegularGridInterpolator((chla_fq, wl_fq), q0_morel, method=method),
            'sq': RegularGridInterpolator((chla_fq, wl_fq), sq_morel, method=method),
            'fq': RegularGridInterpolator((chla_fq, wl_fq), fq_morel, method=method),
            'sfq': RegularGridInterpolator((chla_fq, wl_fq), sfq_morel, method=method)
        }
    
    chla_clip = np.clip(chla, chla_fq.min(), chla_fq.max())
    wl_clip = np.clip(wl, wl_fq.min(), wl_fq.max())
    
    f0 = _MOREL_INTERPOLATORS['f0']([chla_clip, wl_clip])[0]
    sf = _MOREL_INTERPOLATORS['sf']([chla_clip, wl_clip])[0]
    q0 = _MOREL_INTERPOLATORS['q0']([chla_clip, wl_clip])[0]
    sq = _MOREL_INTERPOLATORS['sq']([chla_clip, wl_clip])[0]
    fq = _MOREL_INTERPOLATORS['fq']([chla_clip, wl_clip])[0]
    sfq = _MOREL_INTERPOLATORS['sfq']([chla_clip, wl_clip])[0]
    
    res = [f0, sf, q0, sq, fq, sfq]
    return res


def get_morel_f_over_q(wls, chla, sol_zen_in_degree, method='linear'):
    """
    Get Morel f/Q values for bidirectional effects.
    
    Parameters:
    -----------
    wls : array
        Wavelengths (nm)
    chla : float
        Chlorophyll-a concentration (mg/m3)
    sol_zen_in_degree : float
        Solar zenith angle (degrees)
    method : str
        Interpolation method
        
    Returns:
    --------
    ffqns : array
        f/Q ratios for each wavelength
    """
    sol_zen_in_degree = 30.0
    mu0 = np.cos(sol_zen_in_degree / 180.0 * np.pi)
    
    chla_here = chla if chla >= 0.00001 else 0.1
    
    n_wl = len(wls)
    ffqns = np.zeros(n_wl)
    
    one_minus_mu0 = 1.0 - mu0
    
    for i, wl in enumerate(wls):
        morel_vals = get_morel_numbers(wl, chla_here, method=method)
        [f0, sf, q0, sq, fq, sfq] = morel_vals
        ff = f0 + sf * one_minus_mu0
        qn = q0 + sq * one_minus_mu0
        ffqns[i] = ff / qn
        
    return ffqns


def add_fluorescence(lamb, Rrs, chla, method='semilog'):
    """
    Add chlorophyll fluorescence to Rrs spectrum.
    
    Based on K06 fluorescence analysis with moving peak parameterization.
    Uses Gaussian shape (FWHM = 35 nm) centered at chlorophyll-dependent wavelength.
    
    Parameters:
    -----------
    lamb : array
        Wavelengths (nm)
    Rrs : array
        Base Rrs without fluorescence (sr^-1)
    chla : float
        Chlorophyll-a concentration (mg/m3)
    method : str
        Regression method: 'semilog' (default) or 'linear'
        
    Returns:
    --------
    Rrs_with_fl : array
        Rrs with fluorescence added (sr^-1)
    """
    # Regression coefficients from K06 analysis (eta=0.7, g0=0.0002, nu=0.0)
    # These were derived from 81 field stations
    
    if method == 'semilog':
        # Semi-log regression (recommended)
        # Peak wavelength: wl = p_wl * ln(Chla) + q_wl
        p_wl = 2.807122  # nm per ln(mg/m3)
        q_wl = 682.998020  # nm
        
        # Fluorescence intensity: F = p_fl * ln(Chla) + q_fl
        p_fl = 0.00066110  # sr^-1 per ln(mg/m3)
        q_fl = -0.00111339  # sr^-1
        
        # Calculate peak wavelength and intensity
        if chla > 0:
            wl_peak = p_wl * np.log(chla) + q_wl
            fl_intensity = p_fl * np.log(chla) + q_fl
        else:
            wl_peak = 685.0
            fl_intensity = 0.0
            
    elif method == 'linear':
        # Linear regression (alternative)
        # Peak wavelength: wl = a_wl * Chla + b_wl
        a_wl = 0.11615248  # nm per mg/m3
        b_wl = 686.750118  # nm
        
        # Fluorescence intensity: F = a_fl * Chla + b_fl
        a_fl = 0.0000244659  # sr^-1 per mg/m3
        b_fl = -0.00017439  # sr^-1
        
        wl_peak = a_wl * chla + b_wl
        fl_intensity = a_fl * chla + b_fl
    else:
        raise ValueError(f"Unknown method: {method}. Use 'semilog' or 'linear'")
    
    # Ensure non-negative fluorescence
    if fl_intensity < 0:
        fl_intensity = 0.0
    
    # Gaussian shape (FWHM = 35 nm)
    sigma = 35.0 / 2.355  # Convert FWHM to sigma
    wl_offset = np.arange(-50., 51.)  # -50 to +50 nm range
    gaussian = np.exp(-np.power(wl_offset, 2) / (2 * sigma * sigma))
    
    # Scale by intensity
    fl_spectrum = fl_intensity * gaussian
    
    # Add fluorescence to Rrs
    Rrs_with_fl = Rrs.copy()
    
    # Find center index (closest to peak wavelength)
    idx_center = np.argmin(np.abs(lamb - wl_peak))
    
    # Add fluorescence spectrum
    for i, offset in enumerate(wl_offset):
        idx = idx_center + int(offset)
        if 0 <= idx < len(Rrs_with_fl):
            Rrs_with_fl[idx] += fl_spectrum[i]
    
    return Rrs_with_fl


def get_ccrr_rrs_rt(lamb, chla, mine, cdom, fresnel_add, wl_sky, rho_sky, sol_zen_in_degree,
                    eta=0.7, g0=0.0002, nu=0.0):
    """
    CCRR-RT forward model: Compute remote sensing reflectance (Rrs).
    
    This is the core forward model that computes Rrs spectra from water quality
    parameters using the CCRR-RT (Case 2 Regional CoastColour Radiative Transfer) model.
    
    Parameters:
    -----------
    lamb : array
        Wavelengths (nm)
    chla : float
        Chlorophyll-a concentration (mg/m3)
    mine : float
        Mineral concentration (g/m3)
    cdom : float
        CDOM absorption coefficient at 443nm (1/m)
    fresnel_add : float
        Fresnel reflection correction term
    wl_sky : array
        Wavelengths for rho_sky (nm)
    rho_sky : array
        Sky reflectance spectrum
    sol_zen_in_degree : float
        Solar zenith angle (degrees)
    eta : float, optional
        Pigment backscattering exponent (default: 0.001)
    g0 : float, optional
        Pigment backscattering coefficient (default: 0.001)
    nu : float, optional
        Pigment backscattering spectral slope (default: 0.0)
        
    Returns:
    --------
    Rrs : array
        Remote sensing reflectance (sr^-1)
    """
    global _AB_BOX_CACHE, _AB_INTERP_CACHE, _WAVELENGTH_EXP_CACHE
    
    # Load absorption parameters (cached)
    if _AB_BOX_CACHE is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        dir_AB = os.path.join(current_dir, "independent")
        fname_margalef = "AB_box_Margalef.pkl"
        with open(os.path.join(dir_AB, fname_margalef), 'rb') as f:
            _AB_BOX_CACHE = pickle.load(f)
    
    # Cache AB interpolation results
    lamb_key = _get_lamb_key(lamb)
    if lamb_key not in _AB_INTERP_CACHE:
        wl_ref, AB_box, AB_box_all = _AB_BOX_CACHE
        AA = interpolate_data(wl_ref, AB_box_all[:, 0], lamb)
        BB = interpolate_data(wl_ref, AB_box_all[:, 1], lamb)
        _AB_INTERP_CACHE[lamb_key] = (AA, BB)
    else:
        AA, BB = _AB_INTERP_CACHE[lamb_key]
    
    # Pigment absorption
    as_pig = AA * np.power(chla, -BB)
    a_pig = as_pig * chla
    
    # Cache wavelength-dependent exponential terms
    cache_key = (lamb_key, nu)
    
    if cache_key not in _WAVELENGTH_EXP_CACHE:
        exp_mine = np.exp(-0.01 * (lamb - 440.0))
        exp_cdom = np.exp(-0.017 * (lamb - 400.0))
        lamb_555_pow = np.power(lamb / 555.0, -0.3749)
        lamb_550_pow = np.power(lamb / 550.0, nu)
        
        _WAVELENGTH_EXP_CACHE[cache_key] = (exp_mine, exp_cdom, lamb_555_pow, lamb_550_pow)
    else:
        exp_mine, exp_cdom, lamb_555_pow, lamb_550_pow = _WAVELENGTH_EXP_CACHE[cache_key]
    
    # Mineral absorption (Babin et al., 2003)
    a_mine_440 = 0.025 * mine
    a_mine = a_mine_440 * exp_mine
    
    # CDOM absorption
    a_cdom = cdom * exp_cdom
    
    # Water absorption
    a_wat, b_wat = get_water_ab(lamb)
    
    # Total absorption
    a_tot = a_wat + a_pig + a_mine + a_cdom
    
    # Mineral backscattering
    b_mine_555 = 0.51 * mine
    a_mine_555 = a_mine_440 * np.exp(-0.011 * (555.0 - 440.0))
    g_mine_555 = a_mine_555 + b_mine_555
    g_mine = g_mine_555 * lamb_555_pow
    b_mine = g_mine - a_mine
    bb_mine = 0.019 * b_mine
    
    # Pigment backscattering (Morel 2002 parameterization)
    b_pig_550 = np.power(chla, eta)
    b_pig = b_pig_550 * lamb_550_pow
    bb_pig = g0 * b_pig
    
    # Total backscattering
    bb_tot = bb_pig + bb_mine + 0.5 * b_wat
    
    # Compute above-water Rrs using Morel's f/Q
    chla_morel = min(chla, 10.0)
    f_over_q = get_morel_f_over_q(lamb, chla_morel, sol_zen_in_degree)
    Rrs = 0.53 * f_over_q * bb_tot / (a_tot + bb_tot)
    
    # Add sky reflection contribution
    rho_sky_int = interpolate_data(wl_sky, rho_sky, lamb)
    Rrs += fresnel_add * rho_sky_int
    
    return Rrs


def runForward(wavelengths, chla, mineral, aCDOM, 
               eta_param=0.7, g0_param=0.0002, nu_param=0.0,
               fresnel_add=0.02, rho_sky=None, sol_zen=30.0,
               add_fluorescence_flag=True, fluorescence_method='semilog'):
    """
    Run the CCRR-RT forward model to compute Rrs spectra.
    
    This is the main user-facing function. Given water quality parameters,
    it computes the remote sensing reflectance (Rrs) spectrum.
    
    Parameters:
    -----------
    wavelengths : array-like
        Wavelengths for which to compute Rrs (nm)
        Recommended range: 400-750 nm
        
    chla : float
        Chlorophyll-a concentration (mg/m3)
        Typical range: 0.01 - 100 mg/m3
        
    mineral : float
        Mineral (suspended sediment) concentration (g/m3)
        Typical range: 0.1 - 10 g/m3
        
    aCDOM : float
        CDOM absorption coefficient at 443nm (1/m)
        Typical range: 0.001 - 0.5 (1/m)
        
    eta_param : float, optional
        Pigment backscattering exponent (default: 0.001)
        Typical range: 0.0001 - 0.05
        
    g0_param : float, optional
        Pigment backscattering coefficient (default: 0.001)
        Typical range: 0.0001 - 0.01
        
    nu_param : float, optional
        Pigment backscattering spectral slope (default: 0.0)
        Typical range: 0.0 - 0.1
        
    fresnel_add : float, optional
        Fresnel reflection correction (default: 0.02)
        
    rho_sky : array-like or None, optional
        Sky reflectance spectrum. If None, assumed to be zero.
        Should have same length as wavelengths.
        
    sol_zen : float, optional
        Solar zenith angle in degrees (default: 30.0)
        Note: Currently fixed at 60 degrees in the model
        
    add_fluorescence_flag : bool, optional
        Whether to add chlorophyll fluorescence (default: True)
        Fluorescence is parameterized based on K06 analysis
        
    fluorescence_method : str, optional
        Method for fluorescence calculation: 'semilog' or 'linear' (default: 'semilog')
        'semilog': Fluorescence intensity = p * ln(Chla) + q (recommended)
        'linear': Fluorescence intensity = a * Chla + b
        
    Returns:
    --------
    Rrs : ndarray
        Remote sensing reflectance (sr^-1) at each wavelength
        
    Example:
    --------
    >>> import numpy as np
    >>> import MargModel
    >>> 
    >>> # Define wavelength range
    >>> wavelengths = np.arange(400., 750.)
    >>> 
    >>> # Set water quality parameters
    >>> chla = 2.0       # 2 mg/m3 chlorophyll
    >>> mineral = 0.5    # 0.5 g/m3 suspended sediment
    >>> aCDOM = 0.02     # 0.02 (1/m) CDOM absorption at 443nm
    >>> 
    >>> # Compute Rrs
    >>> Rrs = MargModel.runForward(wavelengths, chla, mineral, aCDOM)
    >>> 
    >>> # Plot results
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(wavelengths, Rrs * 1000)  # Convert to 10^-3 sr^-1
    >>> plt.xlabel('Wavelength (nm)')
    >>> plt.ylabel('Rrs (x10^-3 sr^-1)')
    >>> plt.show()
    """
    # Convert wavelengths to numpy array if needed
    lamb = np.asarray(wavelengths, dtype=float)
    
    # Handle rho_sky
    if rho_sky is None:
        # If no sky reflectance provided, use zeros
        wl_sky = lamb.copy()
        rho_sky = np.zeros_like(lamb)
    else:
        wl_sky = lamb.copy()
        rho_sky = np.asarray(rho_sky, dtype=float)
        
        # Check dimensions
        if len(rho_sky) != len(lamb):
            raise ValueError(f"rho_sky length ({len(rho_sky)}) must match wavelengths length ({len(lamb)})")
    
    # Call the forward model with parameters
    Rrs = get_ccrr_rrs_rt(lamb, chla, mineral, aCDOM, 
                          fresnel_add, wl_sky, rho_sky, sol_zen,
                          eta=eta_param, g0=g0_param, nu=nu_param)
    
    # Add fluorescence if requested
    if add_fluorescence_flag:
        Rrs = add_fluorescence(lamb, Rrs, chla, method=fluorescence_method)
    
    return Rrs


def getDefaultParameters():
    """
    Get default model parameters.
    
    Returns:
    --------
    dict : Dictionary containing default eta, g0, nu values
    """
    return {'eta': 0.7, 'g0': 0.0002, 'nu': 0.0}


if __name__ == "__main__":
    """
    Example usage and testing
    """
    import matplotlib.pyplot as plt
    
    print("=" * 60)
    print("MargModel.py - CCRR-RT Forward Model")
    print("=" * 60)
    
    # Define wavelength range
    wavelengths = np.arange(400., 750.)
    
    # Test Case 1: Typical coastal water
    print("\nTest Case 1: Typical coastal water")
    print("-" * 60)
    chla = 2.0
    mineral = 0.5
    aCDOM = 0.02
    
    print(f"Parameters:")
    print(f"  Chla: {chla} mg/m3")
    print(f"  Mineral: {mineral} g/m3")
    print(f"  aCDOM(443): {aCDOM} (1/m)")
    
    Rrs1 = runForward(wavelengths, chla, mineral, aCDOM)
    print(f"  Max Rrs: {np.max(Rrs1)*1000:.3f} x10^-3 sr^-1")
    print(f"  Rrs at 555nm: {Rrs1[155]*1000:.3f} x10^-3 sr^-1")
    
    # Test Case 2: Clear water (low Chla)
    print("\nTest Case 2: Clear water")
    print("-" * 60)
    chla = 0.5
    mineral = 0.1
    aCDOM = 0.005
    
    print(f"Parameters:")
    print(f"  Chla: {chla} mg/m3")
    print(f"  Mineral: {mineral} g/m3")
    print(f"  aCDOM(443): {aCDOM} (1/m)")
    
    Rrs2 = runForward(wavelengths, chla, mineral, aCDOM)
    print(f"  Max Rrs: {np.max(Rrs2)*1000:.3f} x10^-3 sr^-1")
    print(f"  Rrs at 555nm: {Rrs2[155]*1000:.3f} x10^-3 sr^-1")
    
    # Test Case 3: Turbid water (high mineral)
    print("\nTest Case 3: Turbid water")
    print("-" * 60)
    chla = 5.0
    mineral = 2.0
    aCDOM = 0.05
    
    print(f"Parameters:")
    print(f"  Chla: {chla} mg/m3")
    print(f"  Mineral: {mineral} g/m3")
    print(f"  aCDOM(443): {aCDOM} (1/m)")
    
    Rrs3 = runForward(wavelengths, chla, mineral, aCDOM)
    print(f"  Max Rrs: {np.max(Rrs3)*1000:.3f} x10^-3 sr^-1")
    print(f"  Rrs at 555nm: {Rrs3[155]*1000:.3f} x10^-3 sr^-1")
    
    # Plot comparison
    print("\n" + "=" * 60)
    print("Generating comparison plots...")
    print("=" * 60)
    
    # Figure 1: Original 3 test cases
    plt.figure(figsize=(10, 6))
    plt.plot(wavelengths, Rrs1 * 1000, 'b-', linewidth=2, label='Case 1: Typical coastal')
    plt.plot(wavelengths, Rrs2 * 1000, 'g-', linewidth=2, label='Case 2: Clear water')
    plt.plot(wavelengths, Rrs3 * 1000, 'r-', linewidth=2, label='Case 3: Turbid water')
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Rrs (x10$^{-3}$ sr$^{-1}$)', fontsize=12)
    plt.title('CCRR-RT Forward Model: Example Rrs Spectra', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('MargModel_test_output.png', dpi=150, bbox_inches='tight')
    print("\nPlot 1 saved as: MargModel_test_output.png")
    
    # Figure 2: High Chlorophyll sensitivity (Chla = 10, 50, 100 mg/m3)
    print("\nTest Case 4: High Chlorophyll Sensitivity (with Fluorescence)")
    print("-" * 60)
    mineral_fixed = 0.5
    aCDOM_fixed = 0.01
    chla_values = [10.0, 50.0, 100.0]
    
    print(f"Fixed parameters:")
    print(f"  Mineral: {mineral_fixed} g/m3")
    print(f"  aCDOM(443): {aCDOM_fixed} (1/m)")
    print(f"\nVarying Chlorophyll-a (with fluorescence):")
    
    plt.figure(figsize=(10, 6))
    colors = ['#2ca02c', '#ff7f0e', '#d62728']
    
    for i, chla_val in enumerate(chla_values):
        Rrs_high = runForward(wavelengths, chla_val, mineral_fixed, aCDOM_fixed)
        plt.plot(wavelengths, Rrs_high * 1000, linewidth=2.5, 
                color=colors[i], label=f'Chla = {chla_val:.0f} mg/m³')
        
        max_rrs = np.max(Rrs_high) * 1000
        max_wl = wavelengths[np.argmax(Rrs_high)]
        
        # Check for fluorescence peak around 685 nm
        idx_fl = np.argmin(np.abs(wavelengths - 685))
        rrs_fl = Rrs_high[idx_fl] * 1000
        print(f"  Chla = {chla_val:5.0f} mg/m3: Max Rrs = {max_rrs:.3f} x10^-3 sr^-1 at {max_wl:.0f} nm, Rrs(685nm) = {rrs_fl:.3f}")
    
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Rrs (x10$^{-3}$ sr$^{-1}$)', fontsize=12)
    plt.title('High Chlorophyll Sensitivity with Fluorescence\n(TSM=0.5 g/m³, CDOM=0.01 1/m)', 
             fontsize=14, fontweight='bold')
    plt.legend(fontsize=11, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('MargModel_high_chla_sensitivity.png', dpi=150, bbox_inches='tight')
    print("\nPlot 2 saved as: MargModel_high_chla_sensitivity.png")
    
    # Figure 3: Fluorescence comparison (with vs without)
    print("\nTest Case 5: Fluorescence Effect Comparison")
    print("-" * 60)
    chla_demo = 50.0
    
    print(f"Demonstrating fluorescence effect for Chla = {chla_demo} mg/m3")
    
    Rrs_no_fl = runForward(wavelengths, chla_demo, mineral_fixed, aCDOM_fixed, 
                           add_fluorescence_flag=False)
    Rrs_with_fl = runForward(wavelengths, chla_demo, mineral_fixed, aCDOM_fixed, 
                             add_fluorescence_flag=True)
    
    plt.figure(figsize=(10, 6))
    plt.plot(wavelengths, Rrs_no_fl * 1000, 'b-', linewidth=2.5, label='Without Fluorescence')
    plt.plot(wavelengths, Rrs_with_fl * 1000, 'r-', linewidth=2.5, label='With Fluorescence')
    plt.plot(wavelengths, (Rrs_with_fl - Rrs_no_fl) * 1000, 'g--', linewidth=2, 
             label='Fluorescence Signal')
    
    # Mark fluorescence peak
    fl_peak_wl = 3.088296 * np.log(chla_demo) + 682.024811
    plt.axvline(fl_peak_wl, color='gray', linestyle=':', linewidth=1.5, 
                label=f'FL Peak ({fl_peak_wl:.1f} nm)')
    
    plt.xlabel('Wavelength (nm)', fontsize=12)
    plt.ylabel('Rrs (x10$^{-3}$ sr$^{-1}$)', fontsize=12)
    plt.title(f'Fluorescence Effect (Chla = {chla_demo} mg/m³)', 
             fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.xlim(400, 750)
    plt.tight_layout()
    plt.savefig('MargModel_fluorescence_effect.png', dpi=150, bbox_inches='tight')
    print("\nPlot 3 saved as: MargModel_fluorescence_effect.png")
    
    # Print fluorescence statistics
    fl_signal = Rrs_with_fl - Rrs_no_fl
    max_fl = np.max(fl_signal) * 1000
    max_fl_wl = wavelengths[np.argmax(fl_signal)]
    print(f"  Fluorescence peak: {max_fl:.4f} x10^-3 sr^-1 at {max_fl_wl:.0f} nm")
    print(f"  Predicted peak wavelength: {fl_peak_wl:.1f} nm")
    
    plt.show()
    
    print("\n" + "=" * 60)
    print("All tests completed successfully!")
    print("=" * 60)

