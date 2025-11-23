% Compute global coherence using the methods described by Cimenser et al.
% (2011), "Tracking brain states under general anesthesia by using global
% coherence analysis" - doi: 10.1073/pnas.1017041108
% Inputs:
% data - sources x time matrix
% fs - sampling rate (Hz)
% winSize - size of windows for computing cross power spectral density
% overlap - amount of overlap between windows
% segmentLength - length of epochs in data
% Outputs:
% GC_epoch - frequencies x epochs matrix containing the global coherence
% (GC). num_freqs is determined by nFFT in CPSD calculation.
% LE_epoch - sources x frequencies x epochs matrix containing the leading
% eigenvector.

function [GC_epoch,LE_epoch] = computeGlobalCoherence_epoch(data, fs, winSize, overlap, segmentLength)

nFFT = winSize; % for cpsd later on
nFreqs = nFFT/2+1;

% reshape Y into a 3D matrix
[p1, Tm]            = size(data);
numFullSegments     = floor(Tm / segmentLength);
validTimepoints     = numFullSegments * segmentLength;
dataTrunc           = data(:, 1:validTimepoints); % Truncate the data to fit whole segments only
data_3D             = permute( reshape(dataTrunc, [p1, segmentLength, numFullSegments]), [1 3 2] );
GC_epoch            = zeros(nFreqs,numFullSegments);
LE_epoch            = zeros(p1,nFreqs,numFullSegments);

for ep = 1:numFullSegments

    data_epoch = squeeze(data_3D(:,ep,:));
    nSources = size(data_epoch, 1);
    window = hamming(winSize);
    noverlap = overlap;
    CSD = zeros(nSources, nSources, nFreqs);
    
    % Compute cross-spectral matrix at each frequency
    for i = 1:nSources
        for j = i:nSources
            [Pxy, ~] = cpsd(data_epoch(i,:), data_epoch(j,:), window, noverlap, nFFT, fs);
            CSD(i,j,:) = Pxy;
            if i ~= j
                CSD(j,i,:) = conj(Pxy);  % Ensure Hermitian symmetry
            end
        end
    end
    
    for f = 1:nFreqs
        C = squeeze(CSD(:,:,f));
        [U,D] = eig(C);
        [eigVals,eigVals_idx] = sort(real(diag(D)), 'descend');
        GC_epoch(f,ep) = eigVals(1) / sum(eigVals);
        LE_epoch(:,f,ep) = U(:,eigVals_idx(1));
    end
end