function DAMAS_result = DAMAS(DAS_result, hn, maxIter)
%
% This code implements the DAMAS algorithm
%
% More information about DAMAS can be found in the paper:
%    Brooks, Thomas F and Humphreys, William M, 
%    "A deconvolution approach for the mapping of acoustic sources (DAMAS) determined from phased microphone arrays", 
%    Journal of sound and vibration, 2006.
%
%
% Inputs:
%    DAS_result:  beamforming map, obtained by DAS
%    hn:  steering vector
%    maxIter: the maximum allowable iterations
%
% Outputs:
%    DAMAS_result:  beamforming map, obtained by DAMAS
%
% Author: Hao Liang 
% Last modified by: 21/09/07
%

temp = real(DAS_result);
deps = 0.1;

% Form the dictionary A
[N_X, N_Y, N_mic] = size(hn);
en = reshape(hn, N_X*N_Y, N_mic);
A = (abs(conj(en)*(1./en)').^2)./N_mic^2;
% This is different from the definition of A in https://github.com/jorgengrythe/beamforming (as follow), 
% and I think it may have some minor mistakes... 
% A = (abs(en*en').^2)./N_mic^2; 

% Initialization
Q = zeros(size(temp)); Q0 = temp;

% Solve the inverse problem Y = AQ for Q using Gauss-Seidel iteration 
% -- Y: DAS result
% -- A: reconstruction dictionary
% -- Q: DAMAS_result
for i = 1:maxIter
    for n = 1:N_X*N_Y
        Q(n) = max(0, temp(n) - A(n, 1:n-1)*Q(1:n-1)' ...
            - A(n, n+1:end)*Q0(n+1:end)');
    end
    
    dX = (Q - Q0);
    maxd = max(abs(dX(:)))/mean(Q0(:));
    
    if  maxd < deps
        break;
    end
    
    Q0 = Q;
end

DAMAS_result = Q;

end
