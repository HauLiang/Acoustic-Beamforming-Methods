function CLEAN_SC_result = CLEAN_SC(loopgain, maxIter, CSM, hn)
%
% This code implements the CLEAN-SC algorithm
%
% More information about CLEAN-SC can be found in the paper:
%    Sijtsma, Pieter, 
%    "CLEAN based on spatial source coherence", 
%    International journal of aeroacoustics, 2007.
%
%
% Inputs:
%    loopgain:  loop gain
%    maxIter:   the maximum allowable iterations
%    CSM:  cross-spectrum matrix (CSM)
%    hn:   steering vector
%    z0:   source distance to array
%    N:    number of grid points in each dim 
%
% Outputs:
%    CLEAN_SC_result:  beamforming map, obtained by CLEAN-SC
%
% Author: Hao Liang 
% Last modified by: 21/09/13
%

% Straighten the steering vector
[N_X, N_Y, N_mic] = size(hn);
h = zeros(N_mic, N_X*N_Y);
for i = 1:N_mic
    hn(:,:,i) = flip(hn(:,:,i),1);
    hn(:,:,i) = flip(hn(:,:,i),2);
    h(i,:) = reshape(hn(:,:,i).',1,[])./N_mic;
end

% Scan points setting         
N_scanpoints = N_X*N_Y;

% Start CLEAN SC procedure
Clean_map = zeros(1, N_scanpoints); Degraded_CSM = conj(CSM); Dirty_map = sum(conj(h).*(Degraded_CSM*h), 1);
Dcurr = sum(abs(Degraded_CSM(:))); count = 0; Dprev = 1e8;

while ( Dcurr < Dprev ) && (count < maxIter)
    
    % Determine peak source
    [Map_max, index_max] = max(abs(Dirty_map)); 
    Map_maxCleanBeam = zeros(size(Clean_map)); Map_maxCleanBeam(index_max) = 1;
    ispositiv = real(Dirty_map(index_max)) > 0;  % Delete negative pressure maps
    gmax = Degraded_CSM*h(:,index_max)/Map_max;
    
    % Update degraded CSM 
    Cource = Map_max*(gmax)*gmax';
    Degraded_CSM = Degraded_CSM - loopgain*Cource;
    
    % Update dirty map 
    Dirty_map = sum(conj(h).*(Degraded_CSM*h), 1);

    % Update clean map 
    Clean_map = Clean_map + loopgain*Map_max*Map_maxCleanBeam*ispositiv;

    Dprev = Dcurr; Dcurr = sum(abs(Degraded_CSM(:)));
    count = count + 1;

end


% Delete negative power
Dirty_map(real(Dirty_map)<0) = 0;

% Final beamforming map equal to the sum of clean map and dirty map
CLEAN_SC_result = Clean_map + Dirty_map;
CLEAN_SC_result = reshape(CLEAN_SC_result, N_X, N_Y).';

end


% Delete negative power
Dirty_map(real(Dirty_map)<0) = 0;

% Final beamforming map equal to the sum of clean map and dirty map
CLEAN_SC_result = Clean_map + Dirty_map;
CLEAN_SC_result = reshape(CLEAN_SC_result, N_X, N_Y).';

end
