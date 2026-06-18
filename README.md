# SSOP Multi-Frequency Depth-Sensitive Imaging

**Stony Brook University | NIH R01-Funded | bioRxiv 2026**

Single-Snapshot Optical Properties (SSOP) via structured illumination and Fourier-domain demodulation — enabling depth-resolved tissue characterization from a single image per frequency.

## Overview

Unlike standard 3-phase SFDI, SSOP extracts both AC and DC reflectance from a **single snapshot** using 2D Fourier filtering. This halves acquisition time and enables higher temporal resolution for dynamic in vivo imaging.

Multi-frequency acquisition (5 spatial frequencies) provides depth-sensitivity: high frequencies probe shallow tissue layers, low frequencies probe deeper layers.

## Algorithm

```
Single structured-light image
     ↓  2D FFT
     ↓  Gaussian bandpass filter centered at ±fx  →  AC carrier
     ↓  Low-pass filter  →  DC component
     ↓  Inverse FFT + demodulation → |I_AC|, I_DC
     ↓  Calibration against known phantom
     ↓  LUT inversion → μa(x,y), μs'(x,y) per frequency
```

## Files

| File | Language | Description |
|------|----------|-------------|
| `ssop_demodulation.m` | MATLAB | Fourier-domain SSOP demodulation + multi-frequency pipeline |

## Quick Start

```matlab
% Run the built-in demo
run_ssop_demo()
```

## Depth Sensitivity

Approximate effective probing depth:
```
δ(fx) ≈ 1 / (2 · √(μeff² + (2π·fx)²))
```

Higher spatial frequency → shallower probing depth → layer-resolved tissue imaging.

## Publications

- **Ahmmed et al.** *Depth-Sensitive Optical Property Characterization Using Multi-Frequency Laparoscopic SFDI*, bioRxiv 2026. https://doi.org/10.64898/2026.02.04.703750
- Vervandier & Gioux, Biomed Opt Express 4(10):2040 (2013) – original SSOP

## Author

Rasel Ahmmed | rasel.ahmmed@stonybrook.edu | [Portfolio](https://raselece25.github.io)
