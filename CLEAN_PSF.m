function CLEAN_PSF_result = CLEAN_PSF(loopgain, maxIter, CSM, hn)
%
% This code implements the CLEAN-PSF algorithm
%
% More information about CLEAN-PSF can be found in the paper:
%    Högbom, JA, 
%    "Aperture synthesis with a non-regular distribution of interferometer baselines", 
%    Astronomy and Astrophysics Supplement Series, 1974.
%
%
% Inputs:
%    loopgain:  loop gain
%    maxIter:   the maximum allowable iterations
%    CSM:  cross-spectrum matrix (CSM)
%    hn:   weighted steering vector
%
% Outputs:
%    CLEAN_PSF_result:  beamforming map, obtained by CLEAN-PSF
%
% Author: Hao Liang 
% Last modified by: 21/10/31
%


% Straighten the steering vector
[N_X, N_Y, N_mic] = size(hn); 
h = reshape(hn,[], N_mic).'./N_mic;

% Scan points setting 
N_scanpoints = N_X*N_Y;

% Start CLEAN-PSF procedure
Clean_map = zeros(1, N_scanpoints); Degraded_CSM = CSM; Dirty_map = sum(conj(h).*(Degraded_CSM*h), 1);
Dcurr = sum(abs(Degraded_CSM(:))); count = 0; Dprev = 1e8;
 
while ( Dprev > Dcurr ) && (count < maxIter)
    
    % Determine peak source
    [Map_max, index_max] = max(abs(Dirty_map)); 
    PmaxCleanBeam = zeros(size(Clean_map)); PmaxCleanBeam(index_max) = 1;
    ispositiv = real(Dirty_map(index_max)) > 0;  % Delete negative pressure maps
    
    % Steering vector according to maximum point
    gmax = conj(1./h(:, index_max))*N_mic;
    
    % Update degraded CSM
    Degraded_CSM = Degraded_CSM - loopgain*Map_max*(gmax*gmax');
    
    % Update dirty map 
    Dirty_map = sum(conj(h).*(Degraded_CSM*h), 1);
    
    % Update clean map
    Clean_map = Clean_map + loopgain*Map_max*PmaxCleanBeam*ispositiv;
    
    % Stopping criteria
    Dprev = Dcurr; Dcurr = sum(abs(Degraded_CSM(:)));
    count = count + 1;
    
end

% Delete negative power
Dirty_map(real(Dirty_map)<0) = 0;

% Final beamforming map equal to the sum of clean map and dirty map
CLEAN_PSF_result = Clean_map + Dirty_map;
CLEAN_PSF_result = reshape(CLEAN_PSF_result, N_X, N_Y);

end
