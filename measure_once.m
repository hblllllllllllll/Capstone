%% Dual-Channel ADC: ch0 = mixer_LPF_out, ch1 = Vtune
% Connection: mixer_LPF_out → AI0,  Vtune → AI2,  GND → AGND
function error_value = measure_once()
%% ===== Parameters (modifiable) =====
boardID        = "Board0";
ai_ch_mix      = "ai0";      % ch0: Mixer + LPF output
ai_ch_vtune    = "ai2";      % ch1: Vtune (used for segmentation)
Fs             = 3.33e5;        % Sampling rate (two channels: ≤ 500 kS/s)
captureTime    = 0.0020;     % Acquisition duration (s), should cover multiple sweeps

% ==== NEW: Cancel smoothing and change to extreme value segmentation + head and tail clipping ====
min_sweep_len  = 2;         % Minimum segment length (samples) to avoid false triggers (half-sweep)
pad_mult       = 4;          % Zero-padding multiplier for FFT
n_peaks        = 6;          % Number of peaks to extract per segment
trim_mode      = 'frac';     % 'frac' or 'samples'
trim_head      = 0.15;       % cut head
trim_tail      = 0.12;       % cut tail

% show one frequency
seg_idx_to_show = [];        % 1,2,3,...choose one

%% ===== FMCW params for distance axis =====
VCO_Max_Freq    = 2.25e9;   % [Hz]
VCO_Min_Freq    = 2.09e9;   % [Hz]
Vtune_Freq      = 9e3;      % [Hz]
Vtune_DutyCycle = 0.9;      % duty of up-sweep
VF              = 0.65;     % cable velocity factor
c               = 3e8;      % [m/s]

B    = VCO_Max_Freq - VCO_Min_Freq;      % [Hz] sweep bandwidth
T    = 1 / Vtune_Freq;                    % [s]  period
T_up = T * Vtune_DutyCycle;               % [s]  up-sweep duration
S    = B / T_up;                          % [Hz/s] sweep slope
c_eff = VF * c;                           % [m/s] effective wave speed

%% ===== DAQ Configuration =====
dq = daq("mcc");
addinput(dq, boardID, ai_ch_mix,   "Voltage");  % AI0
addinput(dq, boardID, ai_ch_vtune, "Voltage");  % AI2
dq.Rate = Fs;

%% ===== Data Acquisition =====
data = read(dq, seconds(captureTime));
t    = seconds(data.Time);
y0   = data{:,1};   % ch0: mixer_LPF_out
y1   = data{:,2};   % ch1: Vtune

% Remove DC offset
y0 = y0 - mean(y0);
y1 = y1 - mean(y1);

%% ===== Digital Bandstop Filter: 9 kHz Narrowband Notch (IIR + Zero Phase)=====
f0    = 7e3;      % Notch center frequency (Hz)
bwHz  = 6000;      % Notch bandwidth (Hz)
W0    = f0 / (Fs/2);       % Normalized center frequency (relative to Nyquist)
BW    = bwHz / (Fs/2);     % Normalized bandwidth
[b,a] = iirnotch(W0, BW);  % Second-order notch filter

% Zero-phase bidirectional filtering: avoids phase distortion
y0_f = filtfilt(b, a, y0);
f0    = 10e3;      % Notch center frequency (Hz)
bwHz  = 6000;      % Notch bandwidth (Hz)
W0    = f0 / (Fs/2);       % Normalized center frequency (relative to Nyquist)
BW    = bwHz / (Fs/2);     % Normalized bandwidth
[b,a] = iirnotch(W0, BW);  % Second-order notch filter
y0_f = filtfilt(b, a, y0_f);

f0    = 4e3;      % Notch center frequency (Hz)
bwHz  = 1000;      % Notch bandwidth (Hz)
W0    = f0 / (Fs/2);       % Normalized center frequency (relative to Nyquist)
BW    = bwHz / (Fs/2);     % Normalized bandwidth
[b,a] = iirnotch(W0, BW);  % Second-order notch filter
y0_f = filtfilt(b, a, y0_f);

f0    = 6e3;      % Notch center frequency (Hz)
bwHz  = 6000;      % Notch bandwidth (Hz)
W0    = f0 / (Fs/2);       % Normalized center frequency (relative to Nyquist)
BW    = bwHz / (Fs/2);     % Normalized bandwidth
[b,a] = iirnotch(W0, BW);  % Second-order notch filter
y0_f = filtfilt(b, a, y0_f);

% ===== High-pass parameters =====
fc   = 21e3;           % Cutoff frequency (Hz): Filter out the frequencies on its left
ord  = 6;             % Filter order: The higher the order, the steeper the transition (2 to 8 is reasonable)
Wc   = fc/(Fs/2);     % Normalization

[b,a] = butter(ord, Wc, 'high');
y0_f = filtfilt(b,a,y0_f);   % Zero phase to avoid phase distortion

fc   = 5e3;           % Cutoff frequency (Hz): Filter out the frequencies on its left
ord  = 8;             % Filter order: The higher the order, the steeper the transition (2 to 8 is reasonable)
Wc   = fc/(Fs/2);     % Normalization

[b,a] = butter(ord, Wc, 'high');
y0_f = filtfilt(b,a,y0_f);   % Zero phase to avoid phase loss

%% ===== Find segment boundaries using raw Vtune extrema (both rising & falling) =====
[~, locs_max] = findpeaks( y1, 'MinPeakDistance', min_sweep_len);
[~, locs_min] = findpeaks(-y1, 'MinPeakDistance', min_sweep_len);
locs = sort([locs_max; locs_min]);

N = numel(y1);
locs = locs(locs > 1 & locs < N);

seg_starts = [];
seg_ends   = [];
seg_type   = strings(0);  % "up" / "down"

for k = 1:numel(locs)-1
    s = locs(k);
    e = locs(k+1) - 1;
    if e - s + 1 >= min_sweep_len
        seg_starts(end+1) = s;
        seg_ends(end+1)   = e;
        if y1(locs(k+1)) > y1(locs(k))
            seg_type(end+1) = "up";    % rising
        else
            seg_type(end+1) = "down";  % falling
        end
    end
end

%% only up
up_mask    = (seg_type == "up");
seg_starts = seg_starts(up_mask);
seg_ends   = seg_ends(up_mask);
seg_type   = seg_type(up_mask);
numSweeps  = numel(seg_starts);

%% ===== Visualization: time domain + UP boundaries only =====
figure;
subplot(2,1,1);
plot(t, y1); grid on;
title('Vtune (raw, extrema boundaries, UP-only)'); xlabel('Time (s)'); ylabel('V');
hold on;
for k = 1:numSweeps
    xline(t(seg_starts(k)), '--');
end
hold off;

subplot(2,1,2);
plot(t, y0_f); grid on;
title('Mixer LPF Output (ch0) — UP-only segments'); xlabel('Time (s)'); ylabel('V');
hold on;
for k = 1:numSweeps
    xline(t(seg_starts(k)), '--');
end
hold off;

%% ===== Per-segment processing: trimming + windowing + FFT + peak detection (UP-only) =====
peakFreqs_perSweep  = cell(numSweeps,1);
peakAmpsdB_perSweep = cell(numSweeps,1);

if numSweeps > 0
    % Estimated maximum segment length after clipping
    trimmed_lengths = zeros(numSweeps,1);
    for k = 1:numSweeps
        Lraw = seg_ends(k) - seg_starts(k) + 1;
        switch lower(trim_mode)
            case 'frac'
                cutH = floor(Lraw * trim_head);
                cutT = floor(Lraw * trim_tail);
            case 'samples'
                cutH = round(trim_head);
                cutT = round(trim_tail);
            otherwise
                error('trim_mode should be ''frac'' 或 ''samples''');
        end
        a = seg_starts(k) + max(0,cutH);
        b = seg_ends(k)   - max(0,cutT);
        trimmed_lengths(k) = max(0, b - a + 1);
    end
    maxSegLen = max(max(trimmed_lengths), min_sweep_len);
    Npad = 2^nextpow2(maxSegLen * pad_mult);

    % Frequency axis (only positive frequency, excluding DC)
    freq = Fs*(0:Npad-1)/Npad;
    halfMask = 2:floor(Npad/2);
    freqHalf = freq(halfMask);


    % ★ New: Record Xk and Mag of each section (half spectrum, remove DC)
    XkHalfMatrix  = complex(zeros(numSweeps, numel(freqHalf)));
    MagHalfMatrixLin  = zeros(numSweeps, numel(freqHalf));   % ★ Pre-assigned linear amplitude
    MagHalfMatrix = zeros(numSweeps, numel(freqHalf));
    specMatrix = zeros(numSweeps, numel(freqHalf));

    for k = 1:numSweeps
        s_idx = seg_starts(k);
        e_idx = seg_ends(k);

        % === Cut the beginning and end of the segment ===
        Lraw = e_idx - s_idx + 1;
        switch lower(trim_mode)
            case 'frac'
                cutH = floor(Lraw * trim_head);
                cutT = floor(Lraw * trim_tail);
            case 'samples'
                cutH = round(trim_head);
                cutT = round(trim_tail);
        end
        a = s_idx + max(0,cutH);
        b = e_idx - max(0,cutT);

        % If over-cutting results in insufficient length, try loosening it to meet the minimum length
        if b - a + 1 < min_sweep_len
            spare = floor((cutH+cutT)/2);
            a = max(s_idx, a - spare);
            b = min(e_idx, b + spare);
        end
        if b <= a
            specMatrix(k,:) = -inf;
            peakFreqs_perSweep{k}  = [];
            peakAmpsdB_perSweep{k} = [];
            continue;
        end

        xk = y0_f(a:b);
        xk = xk - mean(xk);             % Remove DC within the segment
        
        % *******changeable: hamming(Nk, 'periodic'), blackman(Nk, 'periodic'), 
        % beta = 5;
        % No window: win = ones(Nk, 1);
        Nk = length(xk);
        win = hann(Nk,'periodic');   
        % win = hamming(Nk,'periodic');
        % win = blackman(Nk,'periodic');
        % win = ones(Nk, 1);
        % win = kaiser(Nk, beta);
        xw = xk .* win;

        % Zero-padding / truncate to Npad
        if Nk < Npad
            xw_pad = [xw; zeros(Npad - Nk, 1)];
        else
            xw_pad = xw(1:Npad);
        end

        % FFT
        Xk  = fft(xw_pad);
        Mag = 20*log10(abs(Xk) + eps);

        % Save the magnitude spectrum (positive frequency, DC removed)
        % ★ Newly added: Half spectrum (remove DC)
        Xk_half = Xk(halfMask).';
        Mag_half_lin = abs(Xk_half);
        Mag_half = 20*log10(abs(Xk_half) + eps);

        % Save (keep consistent with the original specMatrix)
        XkHalfMatrix(k,:)  = Xk_half;
        MagHalfMatrixLin(k,:)  = Mag_half_lin;   % ★ New: Linear Amplitude
        MagHalfMatrix(k,:) = Mag_half;
        specMatrix(k,:)    = Mag_half;
        specMatrix(k,:) = Mag(halfMask).';

        % Peak extraction
        [pk, idx] = findpeaks(specMatrix(k,:), 'SortStr','descend', 'NPeaks', n_peaks);
        peakFreqs_perSweep{k}  = freqHalf(idx);
        peakAmpsdB_perSweep{k} = pk;
    end
else
    % Placeholder when there is no rising edge segment
    Npad = 2^nextpow2(min_sweep_len * pad_mult);
    freq = Fs*(0:Npad-1)/Npad;
    halfMask = 2:floor(Npad/2);
    freqHalf = freq(halfMask);
    specMatrix = zeros(0, numel(freqHalf));
end

%% ===== Print results (main peak per UP segment) =====
fprintf('\n===== Dual-Channel Spectrum Results (UP-only, head/tail trimmed) =====\n');
fprintf('Sampling Rate Fs = %.0f Hz, Duration = %.4f s, UP Segments (valid) = %d\n', Fs, captureTime, numSweeps);
for k = 1:numSweeps
    if ~isempty(peakFreqs_perSweep{k})
        f_main = peakFreqs_perSweep{k}(1);
        a_main = peakAmpsdB_perSweep{k}(1);
        fprintf('UP Segment %d: main peak ~ %.2f kHz, %.1f dB\n', k, f_main/1e3, a_main);
    else
        fprintf('UP Segment %d: no significant peak.\n', k);
    end
end

%% ===== Waterfall plot (spectrum per UP segment) =====
figure;
imagesc(freqHalf/1e3, 1:numSweeps, specMatrix);
axis xy; 
cb = colorbar;             
cb.Label.String = 'Amplitude (dB)';  
cb.Label.FontSize = 12;    
cb.Ticks = -50:5:10;
grid on;
xlabel('Frequency (kHz)'); ylabel('Chirp Index (per sweep)');
title('Frequency–Time Map (Pre-segmented Up-chirp Spectrum)');

%% ===== Waterfall plot (spectrum per UP segment) — X axis as Distance =====
% Map the beat frequency axis to the distance axis: L = (f_b * c_eff) / S
distHalf = (freqHalf * c_eff) / S;   % [m]

figure;
imagesc(distHalf, 1:numSweeps, specMatrix);
axis xy; 
cb = colorbar;             
cb.Label.String = 'Amplitude (dB)';  
cb.Label.FontSize = 12;    
cb.Ticks = -50:5:10;
grid on;
xlabel('Distance (m)'); ylabel('Chirp Index (per sweep)');
title('Range–Time Map (Pre-segmented Up-chirp Spectrum)');

%% ===== Single-segment spectrum (one UP segment, e.g., middle) =====
if numSweeps > 0 && ~isempty(specMatrix)
    valid_rows = find(any(isfinite(specMatrix), 2));
    if isempty(valid_rows)
        warning('There are no valid rising edge segments available, so a single-segment spectrum cannot be output.');
    else
        if isempty(seg_idx_to_show)
            seg_pick = valid_rows(ceil(numel(valid_rows)/2));  % Take the "middle valid UP segment"
        else
            seg_pick = seg_idx_to_show;
            if seg_pick < 1 || seg_pick > numSweeps || ~ismember(seg_pick, valid_rows)
                [~, nearestIdx] = min(abs(valid_rows - min(max(seg_pick,1), numSweeps)));
                seg_pick = valid_rows(nearestIdx);
                warning('seg_idx_to_show Invalid, has been changed to the nearest valid UP segment %d。', seg_pick);
            end
        end

        freq_oneSeg_kHz = freqHalf/1e3;
        mag_oneSeg_dB   = specMatrix(seg_pick, :);

        figure;
        plot(freq_oneSeg_kHz, mag_oneSeg_dB); grid on;
        xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
        title(sprintf('Single-Segment Spectrum: UP segment %d', seg_pick));

        assignin('base', 'seg_pick',        seg_pick);
        assignin('base', 'freq_oneSeg_kHz', freq_oneSeg_kHz);
        assignin('base', 'mag_oneSeg_dB',   mag_oneSeg_dB);
    end
end

%% ===== Compute residual errors using (I - P) from Excel =====
% Configuration：different file/sheet name
ip_excel_file   = 'data6_magnitude.xlsx';              % Save the real part of I-P in Excel
ip_sheet_real   = 'IminusP_real';      % I-P real part sheet

% Read and merge into a real matrix (63×63)
IminusP6 = readmatrix(ip_excel_file, 'Sheet', ip_sheet_real);
assert(isequal(size(IminusP6), [63 63]), 'The IP size must be 63×63');

% Configuration：different file/sheet name
ip_excel_file8   = 'data8_magnitude.xlsx';             
ip_sheet_real8   = 'IminusP_real';      

IminusP8 = readmatrix(ip_excel_file8, 'Sheet', ip_sheet_real8);
assert(isequal(size(IminusP8), [63 63]),'The IP size must be 63×63');

% Configuration：different file/sheet name
ip_excel_file12   = 'data12_magnitude.xlsx';              
ip_sheet_real12   = 'IminusP_real';      % I-P real part sheet

IminusP12 = readmatrix(ip_excel_file12, 'Sheet', ip_sheet_real12);
assert(isequal(size(IminusP12), [63 63]),'The IP size must be 63×63');

% Configuration：different file/sheet name
ip_excel_file14   = 'data14_magnitude.xlsx';       
ip_sheet_real14   = 'IminusP_real';      

IminusP14 = readmatrix(ip_excel_file14, 'Sheet', ip_sheet_real14);
assert(isequal(size(IminusP14), [63 63]),'The IP size must be 63×63');

% Extract 63 frequency points from the FFT half spectrum of each UP segment and form a 63×k real matrix.
if numSweeps == 0
    warning('There is no valid UP segment, so the error cannot be calculated.');
else
    binCount = 63;
    startBin = 1;                               
    totalHalfBins = numel(freqHalf);            
    assert(startBin + binCount - 1 <= totalHalfBins, ...
        'The half spectrum length is not enough to capture 63 bins. Please adjust the startBin or FFT settings.');

    idx63 = startBin : (startBin + binCount - 1);   

    Xk_63xk = complex(zeros(63, numSweeps));
    valid_cnt = 0;
    for kseg = 1:numSweeps
        % Here use the previously saved XkHalfMatrix (the half-spectrum complex number of each segment, row vector)
        x_half = MagHalfMatrixLin(kseg, :);         % 1×(halfBins), real number (amplitude)
        if any(isfinite(x_half))
            valid_cnt = valid_cnt + 1;
            Xk_63xk(:, valid_cnt) = x_half(idx63).';   % Take 63 bins, column vector
        end
    end
    if valid_cnt == 0
        warning('The half spectrum length is not enough to capture 63 bins. Please adjust the startBin or FFT settings.');
    else
        Xk_63xk = Xk_63xk(:, 1:valid_cnt);     % Only keep valid columns

        % Multiplying by (I-P) -> a 63×k error matrix
        E6 = IminusP6 * Xk_63xk;             % IminusP: 63×63（real）
        E8 = IminusP8 * Xk_63xk;  
        E12 = IminusP12 * Xk_63xk;  
        E14 = IminusP14 * Xk_63xk;  
        
        % Calculate the 2-norm by column and take the average
        errs6 = vecnorm(E6, 2, 1);          % 1×k，One ||e||_2 per column
        errs8 = vecnorm(E8, 2, 1); 
        errs12 = vecnorm(E12, 2, 1); 
        errs14 = vecnorm(E14, 2, 1); 

        avg_err6  = mean(errs6);
        avg_err8  = mean(errs8);
        avg_err12 = mean(errs12);
        avg_err14 = mean(errs14);

        fprintf('Average error: e6=%.4f, e8=%.4f, e12=%.4f, e14=%.4f\n', ...
            avg_err6, avg_err8, avg_err12, avg_err14);

    end
end
%% function Output
% Combine the four error averages into a vector output
error_value = [avg_err6, avg_err8, avg_err12, avg_err14];

fprintf('>>> One trial finished, errors = [e6=%.4f, e8=%.4f, e12=%.4f, e14=%.4f]\n', ...
        avg_err6, avg_err8, avg_err12, avg_err14);
end

