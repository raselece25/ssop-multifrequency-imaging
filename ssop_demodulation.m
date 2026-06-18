% ssop_demodulation.m
% Single-Snapshot Optical Properties (SSOP) via Fourier-Domain Demodulation
%
% Extracts AC and DC reflectance components from a SINGLE structured-light
% image using 2D Fourier filtering. This enables real-time optical property
% mapping without the 3-phase acquisition required by standard SFDI.
%
% Algorithm:
%   1. 2D FFT of the raw image
%   2. Band-pass filter centered at (fx, 0) to isolate AC carrier
%   3. Inverse FFT → complex AC image → demodulate amplitude
%   4. Low-pass filter → DC image
%   5. Calibrate and invert (as in sfdi_optical_property_extraction.m)
%
% Reference:
%   Vervandier & Gioux, Biomed Opt Express 4(10):2040 (2013)
%   Ahmmed et al., bioRxiv 2026 – Depth-Sensitive Laparoscopic SFDI/SSOP
%
% Author: Rasel Ahmmed, PhD Candidate, Stony Brook University
% Email : rasel.ahmmed@stonybrook.edu

function [I_AC, I_DC] = ssop_demodulation(I_raw, fx, pixel_size_mm, varargin)
% SSOP_DEMODULATION  Extract AC amplitude and DC from a single SSOP image.
%
%   [I_AC, I_DC] = SSOP_DEMODULATION(I_RAW, FX, PIXEL_SIZE_MM)
%
%   Inputs:
%     I_raw          - Raw single-snapshot image [H x W], float64 or uint16
%     fx             - Structured illumination spatial frequency [mm^-1]
%     pixel_size_mm  - Pixel size in mm (for frequency axis calibration)
%
%   Optional name-value pairs:
%     'BandpassWidth'  - Half-width of Gaussian bandpass filter [mm^-1], default 0.05
%     'LowpassCutoff'  - DC low-pass cutoff [mm^-1], default 0.5*fx
%     'Orientation'    - 'horizontal' | 'vertical' (fringe orientation), default 'horizontal'
%     'ShowFFT'        - true/false – plot intermediate FFT, default false
%
%   Outputs:
%     I_AC  - Demodulated AC amplitude [H x W], double
%     I_DC  - DC component [H x W], double

    p = inputParser;
    addRequired(p, 'I_raw');
    addRequired(p, 'fx');
    addRequired(p, 'pixel_size_mm');
    addParameter(p, 'BandpassWidth',  0.05,         @isnumeric);
    addParameter(p, 'LowpassCutoff',  0.5 * fx,     @isnumeric);
    addParameter(p, 'Orientation',    'horizontal',  @ischar);
    addParameter(p, 'ShowFFT',        false,         @islogical);
    parse(p, I_raw, fx, pixel_size_mm, varargin{:});
    opts = p.Results;

    I_raw = double(I_raw);
    [H, W] = size(I_raw);

    %% Build frequency axes
    % Spatial frequencies in mm^-1
    fx_axis = (-W/2 : W/2-1) / (W * pixel_size_mm);  % [1 x W]
    fy_axis = (-H/2 : H/2-1) / (H * pixel_size_mm);  % [H x 1]
    [FX_grid, FY_grid] = meshgrid(fx_axis, fy_axis);

    %% 2D FFT
    F = fftshift(fft2(I_raw));

    if opts.ShowFFT
        figure('Name','SSOP FFT Spectrum');
        imagesc(fx_axis, fy_axis, log10(abs(F) + 1));
        xlabel('f_x (mm^{-1})'); ylabel('f_y (mm^{-1})');
        title('2D FFT Magnitude (log scale)');
        colormap(hot); colorbar; axis xy;
    end

    %% Bandpass filter to isolate AC carrier at ±fx
    if strcmpi(opts.Orientation, 'horizontal')
        carrier_dist = abs(FX_grid) - opts.fx;    % distance from ±fx ridge
    else
        carrier_dist = abs(FY_grid) - opts.fx;
    end
    sigma_bp = opts.BandpassWidth;
    H_bp     = exp(-carrier_dist.^2 / (2 * sigma_bp^2));  % Gaussian bandpass

    % Keep only positive-frequency half (fx > 0)
    if strcmpi(opts.Orientation, 'horizontal')
        H_bp(FX_grid < 0) = 0;
    else
        H_bp(FY_grid < 0) = 0;
    end

    F_ac = F .* H_bp;

    %% Demodulate: shift carrier to DC and take envelope
    % Shift the spectrum from fx to 0 (demodulation)
    if strcmpi(opts.Orientation, 'horizontal')
        shift_vec = exp(-1j * 2 * pi * opts.fx .* ...
                        repmat((0:W-1), H, 1) * pixel_size_mm);
    else
        shift_vec = exp(-1j * 2 * pi * opts.fx .* ...
                        repmat((0:H-1)', 1, W) * pixel_size_mm);
    end

    ac_complex = ifft2(ifftshift(F_ac)) .* shift_vec * 2;  % ×2 for single sideband
    I_AC       = abs(ac_complex);   % amplitude envelope

    %% Low-pass filter for DC component
    if strcmpi(opts.Orientation, 'horizontal')
        freq_dist_dc = sqrt(FX_grid.^2 + FY_grid.^2);
    else
        freq_dist_dc = sqrt(FX_grid.^2 + FY_grid.^2);
    end
    sigma_lp = opts.LowpassCutoff;
    H_lp     = exp(-freq_dist_dc.^2 / (2 * sigma_lp^2));

    F_dc = F .* H_lp;
    I_DC = real(ifft2(ifftshift(F_dc)));

    % Ensure non-negative
    I_AC = max(I_AC, 0);
    I_DC = max(I_DC, 0);
end


%% ── Multi-Frequency SSOP ─────────────────────────────────────────────────────

function [mu_a_maps, mu_s_prime_maps] = ssop_multifrequency(images_cell, ...
    freq_list, pixel_size_mm, phantom_mu_a, phantom_mu_s_prime, n_tissue, varargin)
% SSOP_MULTIFREQUENCY  Process multiple SSOP images at different spatial frequencies.
%
%   Extracts depth-sensitive optical property maps by processing one snapshot
%   per frequency, then stacking the results. Higher frequencies probe shallower
%   tissue depths; lower frequencies probe deeper tissue.
%
%   Inputs:
%     images_cell      - Cell array of raw images {img_fx1, img_fx2, ...}
%     freq_list        - Corresponding spatial frequencies [mm^-1]
%     pixel_size_mm    - Pixel size in mm
%     phantom_mu_a     - Phantom absorption coefficient [mm^-1]
%     phantom_mu_s_prime - Phantom reduced scattering [mm^-1]
%     n_tissue         - Tissue refractive index (default 1.40)
%
%   Outputs:
%     mu_a_maps       - [H x W x n_freqs] absorption maps
%     mu_s_prime_maps - [H x W x n_freqs] scattering maps

    n_freqs = numel(freq_list);
    fprintf('Processing %d spatial frequencies: %s mm^-1\n', n_freqs, ...
        mat2str(freq_list, 3));

    % Demodulate all images first
    I_AC_all = cell(1, n_freqs);
    I_DC_all = cell(1, n_freqs);

    for k = 1:n_freqs
        fprintf('  Demodulating fx = %.3f mm^-1 ...\n', freq_list(k));
        [I_AC_all{k}, I_DC_all{k}] = ssop_demodulation( ...
            images_cell{k}, freq_list(k), pixel_size_mm, varargin{:});
    end

    [H, W] = size(I_AC_all{1});
    mu_a_maps       = zeros(H, W, n_freqs);
    mu_s_prime_maps = zeros(H, W, n_freqs);

    for k = 1:n_freqs
        % Calibrate against phantom at this frequency
        [R_AC_pred, R_DC_pred] = diffuse_reflectance_model( ...
            phantom_mu_a, phantom_mu_s_prime, freq_list(k), n_tissue);

        cal_AC = R_AC_pred / (mean(I_AC_all{k}(:)) + 1e-12);
        cal_DC = R_DC_pred / (mean(I_DC_all{k}(:)) + 1e-12);

        R_AC = I_AC_all{k} * cal_AC;
        R_DC = I_DC_all{k} * cal_DC;

        % Invert optical properties via LUT
        [mu_a_maps(:,:,k), mu_s_prime_maps(:,:,k)] = ...
            lut_invert(R_AC, R_DC, freq_list(k), n_tissue);

        fprintf('    fx=%.3f mm^-1: mean mu_a=%.4f, mean mu_s''=%.4f mm^-1\n', ...
            freq_list(k), mean(mu_a_maps(:,:,k), 'all'), ...
            mean(mu_s_prime_maps(:,:,k), 'all'));
    end
end


function [mu_a_map, mu_s_prime_map] = lut_invert(R_AC, R_DC, fx, n_tissue)
% LUT_INVERT  Simple LUT inversion for a single (R_AC, R_DC) map.
%   (Reuses diffuse_reflectance_model from sfdi_optical_property_extraction.m)

    N = 150;
    mu_a_vec  = logspace(-4, 0, N);
    mu_sp_vec = logspace(-1, 1.5, N);

    LUT_AC = zeros(N, N);
    LUT_DC = zeros(N, N);
    for i = 1:N
        for j = 1:N
            [LUT_AC(i,j), LUT_DC(i,j)] = ...
                diffuse_reflectance_model(mu_a_vec(i), mu_sp_vec(j), fx, n_tissue);
        end
    end

    [H, W] = size(R_AC);
    mu_a_map       = zeros(H, W);
    mu_s_prime_map = zeros(H, W);

    for px = 1:(H*W)
        dist = (LUT_AC - R_AC(px)).^2 + (LUT_DC - R_DC(px)).^2;
        [~, idx] = min(dist(:));
        [ii, jj] = ind2sub([N, N], idx);
        mu_a_map(px)       = mu_a_vec(ii);
        mu_s_prime_map(px) = mu_sp_vec(jj);
    end

    mu_a_map       = reshape(mu_a_map,       H, W);
    mu_s_prime_map = reshape(mu_s_prime_map, H, W);
end


function [R_d_AC, R_d_DC] = diffuse_reflectance_model(mu_a, mu_s_prime, fx, n)
% (Same physics model as in sfdi_optical_property_extraction.m)
    mu_t_prime = mu_a + mu_s_prime;
    mu_eff     = sqrt(3 * mu_a * mu_t_prime);
    mu_eff_fx  = sqrt(mu_eff^2 + (2*pi*fx).^2);
    A = (1 - r_eff(n)) / (2 * (1 + r_eff(n)));
    R_d_DC = A / (1 + (2 * A * mu_eff)    / (3 * mu_t_prime));
    R_d_AC = A / (1 + (2 * A * mu_eff_fx) / (3 * mu_t_prime));
end

function r = r_eff(n)
    r = -1.440 * n^-2 + 0.710 * n^-1 + 0.668 + 0.0636 * n;
end


%% ── Depth-Sensitivity Estimation ─────────────────────────────────────────────

function probing_depth = sfdi_probing_depth(mu_a, mu_s_prime, fx)
% SFDI_PROBING_DEPTH  Estimate effective probing depth for given optical properties and fx.
%
%   Approximate effective depth δ ≈ 1 / (2 * μ_eff_fx)
%   Higher fx → shallower probing depth.
%
%   Reference: Cuccia et al. 2009, Eq. 6

    mu_t_prime  = mu_a + mu_s_prime;
    mu_eff      = sqrt(3 * mu_a * mu_t_prime);
    mu_eff_fx   = sqrt(mu_eff^2 + (2*pi*fx)^2);
    probing_depth = 1 / (2 * mu_eff_fx);  % in mm
end


%% ── Demo ──────────────────────────────────────────────────────────────────────

function run_ssop_demo()
% RUN_SSOP_DEMO  Demonstrate SSOP demodulation and multi-frequency analysis.

    fprintf('\n--- SSOP Multi-Frequency Demo ---\n');

    [H, W] = deal(128, 128);
    pixel_size_mm = 0.05;    % 50 μm / pixel
    freq_list = [0.05, 0.10, 0.15, 0.20, 0.25];   % mm^-1

    % True tissue optical properties (2-layer heterogeneous phantom)
    mu_a_true_surf = 0.015;   mu_sp_true_surf = 1.20;   % superficial layer
    mu_a_true_deep = 0.008;   mu_sp_true_deep = 0.80;   % deeper layer

    % Calibration phantom
    phantom_mu_a  = 0.010;
    phantom_mu_sp = 1.00;
    n_tissue      = 1.40;

    x = linspace(0, W-1, W) * pixel_size_mm;
    [X, ~] = meshgrid(x, x);

    images_cell = cell(1, numel(freq_list));
    for k = 1:numel(freq_list)
        fx_k = freq_list(k);
        [R_ac, R_dc] = diffuse_reflectance_model(mu_a_true_surf, mu_sp_true_surf, ...
                                                  fx_k, n_tissue);
        % Simulate single-snapshot SSOP image
        img = R_dc + R_ac * cos(2*pi*fx_k*X) + 0.002*randn(H,W);
        images_cell{k} = img;
    end

    % Multi-frequency SSOP
    [mu_a_maps, mu_sp_maps] = ssop_multifrequency(images_cell, freq_list, ...
        pixel_size_mm, phantom_mu_a, phantom_mu_sp, n_tissue);

    % Show probing depths
    fprintf('\nProbing depths (for mean tissue optical properties):\n');
    for k = 1:numel(freq_list)
        d = sfdi_probing_depth(mean(mu_a_maps(:,:,k), 'all'), ...
                               mean(mu_sp_maps(:,:,k), 'all'), freq_list(k));
        fprintf('  fx=%.2f mm^-1: probing depth ≈ %.3f mm\n', freq_list(k), d);
    end

    % Visualization
    figure('Name','SSOP Multi-Frequency Optical Property Maps', ...
           'Position', [50 50 1200 400]);
    n = numel(freq_list);
    for k = 1:n
        subplot(2, n, k);
        imagesc(mu_a_maps(:,:,k)); colorbar; axis image off;
        title(sprintf('\\mu_a: fx=%.2f', freq_list(k)));
        colormap(gca, hot);

        subplot(2, n, n+k);
        imagesc(mu_sp_maps(:,:,k)); colorbar; axis image off;
        title(sprintf('\\mu_s'': fx=%.2f', freq_list(k)));
        colormap(gca, parula);
    end
    sgtitle('SSOP Multi-Frequency Depth-Sensitive Optical Property Maps');
end
