% Compute weighted Phase-Lag Index among all channel pairs. 

% data: 2D matrix of size [N_channels x N_timepoints]. 
% fs: sampling frequency in Hz. 
% fpass: 2-element vector [fLow fHigh] specifying passband.
% WPLI_epoch: N_channels x N_channels x num_epochs matrix of wPLI values.

function WPLI_epoch = compute_wPLI_epoch(data, fpass, fs, segmentLength) 

% Reshape data into a 3D matrix
[p1, Tm]            = size(data);
numFullSegments     = floor(Tm / segmentLength);
validTimepoints     = numFullSegments * segmentLength;
dataTrunc           = data(:, 1:validTimepoints); % Truncate the data to fit whole segments only
data_3D             = permute( reshape(dataTrunc, [p1, segmentLength, numFullSegments]), [1 3 2] );
WPLI_epoch          = zeros(p1, p1, numFullSegments);

for ep = 1:numFullSegments

    % Band-pass filter each channel
    bpOrder = 2*round(fs); 
    filtCoeffs = fir1(bpOrder, [fpass(1)/(fs/2), fpass(2)/(fs/2)]); 
    data_filt = filtfilt(filtCoeffs, 1, squeeze(data_3D(:,ep,:)).').'; 
    
    % Compute analytic signal & phase via Hilbert transform
    analyticSignal = hilbert(data_filt.').'; 
   
    % Compute cross-spectral imaginary component across the epoch
    for i = 1:p1
        for j = i+1:p1
            crossSig = analyticSignal(i,:) .* conj(analyticSignal(j,:));
            imagCross = imag(crossSig); 
    
            % Weighted Phase Lag Index
            num = abs(mean(imagCross, 2));    % | E[ Im(X_{ij}(t)) ] |
            den = mean(abs(imagCross), 2);    % E[ |Im(X_{ij}(t))| ]
    
            wpli_val = 0;
            if den > 0
                wpli_val = num / den;
            end
    
            WPLI_epoch(i,j,ep) = wpli_val;
            WPLI_epoch(j,i,ep) = wpli_val; 
        end
    end
    
end