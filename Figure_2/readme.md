## Calcium Analysis Scripts

### `00_intensity.m`
- Allows manual clicking of ROIs in an image stack  
- Extracts and plots green-channel fluorescence over time  
- Saves ROI data  
- Passes extracted data to an additional analysis function  

### `01_intensity_analysis.m`
- Detects calcium-imaging spikes in each ROI  
- Computes spike amplitudes and ΔF/F₀  
- Measures synchronous vs. individual neuron firing rates  
- Generates spike rasters and normalized traces  
- Performs spectral analysis  
- Saves all figures, results, and statistics (including Excel output)

### `02_fft_neuron.m`
- Removes slow baseline trends from fluorescence traces  
- Applies windowing and normalization  
- Performs zero-padded FFT  
- Identifies strongest spectral peaks  
- Returns each peak’s frequency, period, power, and relative power %

### `03_spike_detection.m`
- Thresholds the signal  
- Groups consecutive above-threshold samples into spikes  
- Returns the index of the maximum point within each detected spike  

### `04_cameron_gcamp.ipynb`
- Normalizes each ROI signal  
- Detects spikes and amplitudes via simple thresholding  
- Computes per-neuron and synchronous spike frequencies  
- Calculates pairwise correlation between active ROIs  
- Saves normalized traces  