function MUSIC_result = MUSIC(CSM, hn, nSources)
%
% This code implements the MUSIC algorithm
%
% More information about MUSIC can be found in the paper:
%    Schmidt, Ralph, 
%    "Multiple emitter location and signal parameter estimation", 
%    IEEE transactions on antennas and propagation, 1986.
%
%
% Inputs:
%    CSM:  cross-spectrum matrix (CSM)
%    hn:   weighted steering vector
%    nSources:   number of sources
%
% Outputs:
%    MUSIC_result:  beamforming map, obtained by MUSIC
%
% Author: Hao Liang 
% Last modified by: 21/09/15
%

% Parameters setting
[N_X, N_Y, N_mic] = size(hn);

% Diagonal reload
CSM = CSM + trace(CSM)/(N_mic^2)*eye(N_mic, N_mic);
CSM = CSM/N_mic;

% Eigenvectors of CSM
[Vec, Val] = eig(CSM);
[~, Seq] = sort(max(Val));

% Noise eigenvectors
Vn = Vec(:,Seq(1:end-nSources));

% Spatial spectrum (pseudo spectral) imaging
MUSIC_result = zeros(N_X, N_Y);
for ii = 1:N_X
    for jj = 1:N_Y
        e = reshape(hn(ii,jj,:), N_mic, 1); 
        MUSIC_result(ii,jj) = 1./(e'*(Vn*Vn')*e);
    end
end

end
