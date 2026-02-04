% Requirements
% - WFDB Toolbox (for rdsamp)
% - bandpass_filter_function (own one) 
% - normalisation_function (own one)
% - approximateEntropy Predictive Maintenance Toolbox in R2018a
% - correlationDimension Predictive Maintenance Toolbox in R2018a
% - DFA_fun Downloaded from File_Exchange
% - Higuchi_FD Downloaded from File_Exchange
% - HurstExponent Downloaded from File_Exchange
% - Katz_FD Downloaded from File_Exchange
% - lyapunovExponent Downloaded from File_Exchange

% Clear command window, workspace, and close figures
clc;                    % clear command window
clear;                  % remove all variables from workspace
close all;              % close all figure windows
tic

% ---------- SETTINGS ----------
% Set the top-level folder containing patient subfolders (PTB database)
rootDir = 'C:\Users\Johan\OneDrive - ucp.pt\Desktop\Privat\Portugal\UCP\Semester_1\Signal_and_Image_Processing\Project\ptb-diagnostic-ecg-database-1.0.0\ptb-diagnostic-ecg-database-1.0.0';  % path to PTB DB
% Define filename for the final Excel output
outFile = 'ECG_Features-Final.xlsx';  % output Excel file

% List directory entries in rootDir
patientFolders = dir(fullfile(rootDir,'*'));  % returns struct array of files + folders
% Keep only entries that are directories
patientFolders = patientFolders([patientFolders.isdir]);  % logical indexing to select only directories
% Remove '.' and '..' pseudo-folders from the list
patientFolders = patientFolders(~ismember({patientFolders.name},{'.','..'}));  % filter out '.' and '..'

% Initialize container for flattened rows (one row per patient)
AllRows = [];           % will be N-patients x N-features_after_stats
% Initialize cell array to store patient folder names as IDs
PatientIDs = {};        % store patient folder names

% Loop over each patient folder
for p = 1:length(patientFolders)                 % iterate through patient folders
    folderName = fullfile(rootDir, patientFolders(p).name);  % full path to current patient folder
    datFiles = dir(fullfile(folderName, '*.dat'));           % find.dat files in folder
    if isempty(datFiles)                             % if no .dat file found
        fprintf('No .dat in %s -> skipped\n', folderName);  % report and skip
        continue;                                    % move to next patient folder
    end

    firstDatFile = datFiles(1).name;                 % take the first .dat filename
    recName = firstDatFile(1:end-4);                 % remove '.dat' extension -> record name

    currentFolder = pwd;                             % save current working directory
    cd(folderName);                                  % change directory to patient folder so rdsamp can find record

    try
        [signal, Fs, ~] = rdsamp(recName, [], 10000); % read 10000 samples using WFDB rdsamp (signal: samples x leads, Fs: sampling freq)
    catch ME
        fprintf('Error loading %s: %s\n', recName, ME.message);  % print error message if read fails
        cd(currentFolder);                            % restore working directory
        continue;                                     % skip this patient
    end

    cd(currentFolder);                                % restore the original working directory after loading

    
    leads = 15;                     
    
    % -----------------------------
    % Applying normalization, Bandpass Butterworth Filter and Notch filter
    % -----------------------------
    normalized_signal = normalisation_function(signal);     % function is provided separately
    filtered_signal = bandpass_filter_function(Fs, normalized_signal, leads); % Filter function is provided separately

    % -----------------------------
    % Create 1-second windows and apply DWT
    % -----------------------------
    seconds = 1;                                    % window duration in seconds
    window_size = Fs * seconds;                     % samples per 1-second window (e.g., Fs*1)
    waveletName = 'sym7';                           % wavelet family name (Symlet 7)
    levels = 3;                                     % decomposition level 3 (A3 + D1 + D2 + D3)

    [numSamples, ~] = size(filtered_signal);                % number of samples (rows) in filtered_signal
    numCompleteWindows = floor(numSamples / window_size);   % how many full 1-second windows fit

    % Frequency ranges per detail/approx band (approximate, for Lyapunov input)
    % D1 ~ Fs/4 .. Fs/2, D2 ~ Fs/8 .. Fs/4, D3 ~ Fs/16 .. Fs/8, A3 ~ 0 .. Fs/16
    freq_D1 = [Fs/4, Fs/2];  % frequency band for detail level 1
    freq_D2 = [Fs/8, Fs/4];  % frequency band for detail level 2
    freq_D3 = [Fs/16, Fs/8]; % frequency band for detail level 3
    freq_A3 = [0, Fs/16];    % frequency band for approximation level 3

    epsilon = 1e-12;         % small constant to avoid log(0) or division-by-zero

    % --- preallocate feature matrices (rows = windows, cols = leads) ---
    approx_entro_A3 = zeros(numCompleteWindows,leads); approx_entro_D1 = zeros(numCompleteWindows,leads);  % ApEn matrices
    approx_entro_D2 = zeros(numCompleteWindows,leads); approx_entro_D3 = zeros(numCompleteWindows,leads);
    corrDim_A3 = zeros(numCompleteWindows,leads); corrDim_D1 = zeros(numCompleteWindows,leads);             % Correlation dimension matrices
    corrDim_D2 = zeros(numCompleteWindows,leads); corrDim_D3 = zeros(numCompleteWindows,leads);
    DFA_A3 = zeros(numCompleteWindows,leads); DFA_D1 = zeros(numCompleteWindows,leads);                      % DFA exponent matrices
    DFA_D2 = zeros(numCompleteWindows,leads); DFA_D3 = zeros(numCompleteWindows,leads);
    energy_A3 = zeros(numCompleteWindows,leads); energy_D1 = zeros(numCompleteWindows,leads);                % Energy matrices
    energy_D2 = zeros(numCompleteWindows,leads); energy_D3 = zeros(numCompleteWindows,leads);
    H_A3 = zeros(numCompleteWindows,leads); H_D1 = zeros(numCompleteWindows,leads);                           % Higuchi FD matrices
    H_D2 = zeros(numCompleteWindows,leads); H_D3 = zeros(numCompleteWindows,leads);
    EH_A3 = zeros(numCompleteWindows,leads); EH_D1 = zeros(numCompleteWindows,leads);                         % Hurst exponent matrices
    EH_D2 = zeros(numCompleteWindows,leads); EH_D3 = zeros(numCompleteWindows,leads);
    katzDim_A3 = zeros(numCompleteWindows,leads); katzDim_D1 = zeros(numCompleteWindows,leads);               % Katz FD matrices
    katzDim_D2 = zeros(numCompleteWindows,leads); katzDim_D3 = zeros(numCompleteWindows,leads);
    logEnergy_A3 = zeros(numCompleteWindows,leads); logEnergy_D3 = zeros(numCompleteWindows,leads);           % log energy matrices
    logEnergy_D2 = zeros(numCompleteWindows,leads); logEnergy_D1 = zeros(numCompleteWindows,leads);
    Elay_A3 = zeros(numCompleteWindows,leads); Elay_D1 = zeros(numCompleteWindows,leads);                     % Lyapunov exponent matrices
    Elay_D2 = zeros(numCompleteWindows,leads); Elay_D3 = zeros(numCompleteWindows,leads);
    ShanEn_A3 = zeros(numCompleteWindows,leads); ShanEn_D1 = zeros(numCompleteWindows,leads);                 % Shannon entropy matrices
    ShanEn_D2 = zeros(numCompleteWindows,leads); ShanEn_D3 = zeros(numCompleteWindows,leads);

    % --- loop over leads and windows ---
    for lead = 1:leads                                   % iterate over channels/leads
        x = filtered_signal(:, lead);                    % extract full signal for this lead
        x_trim = x(1:numCompleteWindows * window_size);  % trim signal to integer number of windows
        windows = reshape(x_trim, window_size, numCompleteWindows);  % reshape into [window_size x numWindows] columns

        for w = 1:numCompleteWindows                % iterate over each 1-second window
            segment = windows(:, w);                % extract 1-second segment (column vector)

            % --- 1) DWT decomposition ---
            [c, l] = wavedec(segment, levels, waveletName);  % perform discrete wavelet decomposition
            % wavedec returns vector c of coefficients and bookkeeping vector l
            % wrcoef reconstructs approximation/detail coefficients at same length as segment
            A3 = wrcoef('a', c, l, waveletName, 3);  % reconstruct approximation at level 3 -> A3 (same length as segment)
            D3 = wrcoef('d', c, l, waveletName, 3);  % reconstruct detail level 3
            D2 = wrcoef('d', c, l, waveletName, 2);  % reconstruct detail level 2
            D1 = wrcoef('d', c, l, waveletName, 1);  % reconstruct detail level 1

            % --- 2) Feature extraction ---

            % ---------- Approximate Entropy (ApEn) ----------
            % ApEn: measures regularity and unpredictability of fluctuations.
            % Formula (high-level): ApEn(m,r) = Phi(m,r) - Phi(m+1,r)
            % where Phi(m,r) = (1/(N-m+1)) * sum_{i=1}^{N-m+1} log(C_i^m(r)),
            % and C_i^m(r) counts the fraction of vectors within tolerance r of vector i.
            approx_entro_A3(w,lead) = approximateEntropy(A3);  % user function: approximateEntropy(x) 
            approx_entro_D1(w,lead) = approximateEntropy(D1);  % compute ApEn for D1
            approx_entro_D2(w,lead) = approximateEntropy(D2);  % compute ApEn for D2
            approx_entro_D3(w,lead) = approximateEntropy(D3);  % compute ApEn for D3

            % ---------- Correlation Dimension (Grassberger-Procaccia D2) ----------
            % High-level formula: C(r) = (1/(N*(N-1))) * sum_{i != j} H(r - ||X_i - X_j||)
            % D2 ≈ d log(C(r)) / d log(r) for small r (slope of log-log plot)
            corrDim_A3(w,lead) = correlationDimension(A3);
            corrDim_D1(w,lead) = correlationDimension(D1);
            corrDim_D2(w,lead) = correlationDimension(D2);
            corrDim_D3(w,lead) = correlationDimension(D3);

            % ---------- Detrended Fluctuation Analysis (DFA) ----------
            % DFA: compute F(s) for multiple box sizes s and fit F(s) ~ s^alpha, where alpha = scaling exponent
            % alpha ~ 0.5 => white noise, alpha~1 => 1/f noise, alpha>1 => trending behaviour
            [alphaA3,~] = DFA_fun(A3,[4 8 16 32 64 128],1); 
            DFA_A3(w,lead) = alphaA3(1);                      % store first alpha (or appropriate scale)
            [alphaD1,~] = DFA_fun(D1,[4 8 16 32 64 128],1);
            DFA_D1(w,lead) = alphaD1(1);
            [alphaD2,~] = DFA_fun(D2,[4 8 16 32 64 128],1);
            DFA_D2(w,lead) = alphaD2(1);
            [alphaD3,~] = DFA_fun(D3,[4 8 16 32 64 128],1);
            DFA_D3(w,lead) = alphaD3(1); 

            % ---------- Energy ----------
            % E = sum_{n=1}^N x[n]^2  (signal energy over the window)
            energy_A3(w,lead) = sum(A3.^2);                   % energy of A3 subband in this window
            energy_D1(w,lead) = sum(D1.^2);                   % energy of D1 subband
            energy_D2(w,lead) = sum(D2.^2);                   % energy of D2 subband
            energy_D3(w,lead) = sum(D3.^2);                   % energy of D3 subband

            % ---------- Higuchi Fractal Dimension (HFD) ----------
            % HFD estimates fractal dimension D by constructing k-time series and L(k) curve,
            % then fitting L(k) ~ k^{-D} (slope of log-log gives D).
            H_A3(w,lead) = Higuchi_FD(A3,3);                   
            H_D1(w,lead) = Higuchi_FD(D1,3);
            H_D2(w,lead) = Higuchi_FD(D2,3);
            H_D3(w,lead) = Higuchi_FD(D3,3);

            % ---------- Hurst Exponent ----------
            % H measures long-term memory; one relation: R/S ~ n^H,
            % where R/S is rescaled range for window size n. H in (0,1).
            EH_A3(w,lead) = HurstExponent(A3);                 
            EH_D1(w,lead) = HurstExponent(D1);
            EH_D2(w,lead) = HurstExponent(D2);
            EH_D3(w,lead) = HurstExponent(D3);

            % ---------- Katz Fractal Dimension ----------
            % Katz FD = log10(n) / (log10(n) + log10(d/L)),
            % where n = number of samples, d = max distance from first point,
            % and L = total length (sum of distances between successive points).
            katzDim_A3(w,lead) = Katz_FD(A3,3);                
            katzDim_D1(w,lead) = Katz_FD(D1,3);
            katzDim_D2(w,lead) = Katz_FD(D2,3);
            katzDim_D3(w,lead) = Katz_FD(D3,3);

            % ---------- Logarithmic Energy ----------
            % logEnergy = sum( log2( x[n]^2 + eps ) ) to avoid log(0)
            logEnergy_A3(w,lead) = sum(log2(A3.^2 + eps));    
            logEnergy_D3(w,lead) = sum(log2(D3.^2 + eps));
            logEnergy_D2(w,lead) = sum(log2(D2.^2 + eps));
            logEnergy_D1(w,lead) = sum(log2(D1.^2 + eps));

            % ---------- Lyapunov Exponent (largest) ----------
            % Largest Lyapunov exponent λ estimates average exponential rate of divergence:
            % λ ≈ (1/T) * Σ log( |δx(t+Δ)| / |δx(t)| )
            Elay_A3(w,lead) = lyapunovExponent(A3,freq_A3(1,2));  
            Elay_D1(w,lead) = lyapunovExponent(D1,freq_D1(1,2));
            Elay_D2(w,lead) = lyapunovExponent(D2,freq_D2(1,2));
            Elay_D3(w,lead) = lyapunovExponent(D3,freq_D3(1,2));

            % ---------- Classical Shannon Entropy (based on sample energy distribution) ----------
            % For samples x[n], compute p_i = x[n]^2 / sum(x^2) then H = -Σ p_i log2 p_i
            % This is a discrete Shannon entropy of the normalized energy distribution.
            x2_A3 = A3.^2;                                     % square samples -> energy per sample
            p_A3 = x2_A3 / sum(x2_A3);                         % normalize to probability distribution p_i
            p_A3(p_A3==0) = eps;                               % avoid log2(0) by replacing zeros with eps
            ShanEn_A3(w,lead) = -sum(p_A3 .* log2(p_A3));      % Shannon entropy (bits)

            x2_D1 = D1.^2;                                     % repeat for D1
            p_D1 = x2_D1 / sum(x2_D1);
            p_D1(p_D1==0) = eps;
            ShanEn_D1(w,lead) = -sum(p_D1 .* log2(p_D1));

            x2_D2 = D2.^2;                                     % repeat for D2
            p_D2 = x2_D2 / sum(x2_D2);
            p_D2(p_D2==0) = eps;
            ShanEn_D2(w,lead) = -sum(p_D2 .* log2(p_D2));

            x2_D3 = D3.^2;                                     % repeat for D3
            p_D3 = x2_D3 / sum(x2_D3);
            p_D3(p_D3==0) = eps;
            ShanEn_D3(w,lead) = -sum(p_D3 .* log2(p_D3));

        end % end windows loop
    end % end leads loop

    % -----------------------------
    % Build feature matrix and compute statistics (compress over time)
    % -----------------------------
    bands  = {'A3','D1','D2','D3'};                      
    groups = {'energy','approx_entro','corrDim','DFA','H','EH','katzDim','logEnergy','Elay','ShanEn'};  % feature groups order

    Feature_Matrix = [];                                % initialize empty matrix to concatenate feature blocks
    for g = 1:numel(groups)                             % iterate feature groups in the exact same order as naming block
        for b = 1:numel(bands)                          % iterate subbands in order A3,D1,D2,D3
            varname = sprintf('%s_%s', groups{g}, bands{b});  % variable name string
            M = eval(varname);                          % retrieve variable by name (M is [numWindows x leads])
            if isequal(size(M), [leads, size(M,1)])     % check orientation; sometimes stored transposed
                M = M.';                                % transpose if needed to ensure rows = windows, cols = leads
            end
            Feature_Matrix = [Feature_Matrix, M];      % horizontally concatenate block -> building columns in consistent order
        end
    end

    % Determine feature matrix size: rows=windows, cols= number of base feature columns (feature×band×lead)
    [numWindows, numCols] = size(Feature_Matrix);      % numCols equals groups*bands*leads
    stats_per_col = zeros(numCols, 6);                 % preallocate for 6 summary statistics per column

    for c = 1:numCols                                  % loop each feature column (feature×lead)
        col = Feature_Matrix(:, c);                    % extract time series for this feature column across windows
        col = col(~isnan(col));                        % remove NaNs that may result from failed computations
        if isempty(col)                                % if no valid values after NaN removal
            stats_per_col(c, :) = NaN(1,6);            % set all six stats to NaN
            continue;                                  % skip to next column
        end
        % Compute the 6 summary statistics used to compress time series:
        % mean  = (1/N) * sum_i x_i
        % std   = sqrt( (1/(N-1)) * sum_i (x_i - mean)^2 )
        % p95   = 95th percentile
        % var   = variance (var(col) uses normalization by N-1)
        % median= middle value
        % kurtosis = measure of tailedness
        stats_per_col(c, :) = [mean(col), std(col), prctile(col,95), var(col), median(col), kurtosis(col)];
    end

    % Flatten stats_per_col: reshape to a single row (1 x (numCols*6))
    stats_row = reshape(stats_per_col.', 1, []);       % transpose then reshape to row-major with all stat blocks
    AllRows = [AllRows; stats_row];                    % append this patient as a row (grow AllRows)
    PatientIDs{end+1} = patientFolders(p).name;        % store patient folder name as ID
    fprintf('Processed patient %s\n', patientFolders(p).name);  % progress message

end % end patient loop

%% -----------------------------
% Final z-score normalization across patients 
% -----------------------------
if isempty(AllRows)                                          % check that at least one patient was processed
    error('No valid patients processed. Nothing to save.');  % error if AllRows empty
end

% Compute column-wise mean (mu) and std (sigma) across patients
mu = mean(AllRows, 1, 'omitnan');                   % 1 x numFeatures: mean of each column across patients
sigma = std(AllRows, 0, 1, 'omitnan');              % 1 x numFeatures: std of each column across patients

% avoid division by zero by replacing zeros with 1 (no scaling)
sigma(sigma == 0) = 1;                              % replace zero std with 1 to keep z-scoring stable
AllRows_z = (AllRows - mu) ./ sigma;                % z-score each feature column across patients: (x - mu)/sigma

% Replace AllRows with z-scored matrix for saving
AllRows = AllRows_z;                                % use z-scored features from here on

% Define the 6 statistical labels for later naming
statLabels = {'mean','std','p95','var','median','kurtosis'};  % labels for the 6 summary stats

% ---------------------------------------------
% Columns pattern: <feature>_<subband>_LeadXX_<statistic>
% ---------------------------------------------
statLabels = {'mean','std','p95','var','median','kurtosis'};  
bands  = {'A3','D1','D2','D3'};                               
groups = {'energy','approx_entro','corrDim','DFA','H','EH','katzDim','logEnergy','Elay','ShanEn'};  % group order

% Number of feature columns before statistics expansion 
numFeatCols = size(AllRows,2) / numel(statLabels);  % total columns divided by 6 stats => base-feature columns

% Build base column names in the exact order as Feature_Matrix (group -> band -> lead)
Column_Names = cell(1, numFeatCols);                 % preallocate cell array for base names
idx = 1;                                             % index counter
for g = 1:numel(groups)                              % iterate each feature group in order
    for b = 1:numel(bands)                           % iterate each band in order
        for lead = 1:leads                           % iterate each lead
            Column_Names{idx} = sprintf('%s_%s_Lead%02d', groups{g}, bands{b}, lead);
            idx = idx + 1;                           % increment index
        end
    end
end


% Expand base names with statistics suffixes to create final column names
colNames = reshape(strcat(repelem(Column_Names,1,numel(statLabels)),'_', repmat(statLabels,1,numel(Column_Names))), 1, []);  
colNames = matlab.lang.makeValidName(colNames);       % ensure valid MATLAB variable names

% ---------------------------------------------
% Build table and save to Excel
% ---------------------------------------------
T = array2table(AllRows, 'VariableNames', colNames);  % convert AllRows to table with variable names
T.patient_id = PatientIDs(:);                         % add patient_id column from PatientIDs cell array
T = movevars(T,'patient_id','Before',1);              % move patient_id to the first column

writetable(T, outFile);                               % write table to Excel file (outFile)
fprintf('\nFinished. Saved %d patients to %s\n', size(AllRows,1), outFile);  % final status message
toc                                                   % Just for fun Time measurement -> 174 min.  