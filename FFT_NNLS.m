function FFT_NNLS_result = FFT_NNLS(DAS_result, PSF, maxIter)
%
% This code implements the FFT-NNLS algorithm
%
% More information about FFT-NNLS can be found in the paper:
%    Ehrenfried, Klaus and Koop, Lars, 
%    "Comparison of iterative deconvolution algorithms for the mapping of acoustic sources",
%    AIAA journal, 2007.
%
%
% Inputs:
%    DAS_result:  beamforming map, obtained by DAS
%    PSF:  point spread function (PSF)
%    maxIter: the maximum allowable iterations
%
% Outputs:
%    FFT_NNLS_result:  beamforming map, obtained by FFT-NNLS
%
% Author: Hao Liang 
% Last modified by: 21/09/08
%

% To avoid wraparound effects, zero padding is used
DAS_zeropad = real(zeropad(DAS_result)); PSF_zeropad = zeropad(PSF);
FFT_NNLS_result_zeropad = zeropad(zeros(size(DAS_result,1)));

% Precompute fft of PSF
fft_PSF = fft2(PSF_zeropad);
   
n = 0;
while n < maxIter
    
    n = n + 1;
    
    % Calculate the residual vector and gradient
    [r, grad] = nnls(DAS_zeropad,FFT_NNLS_result_zeropad,fft_PSF);  
    
    % Compute the vector w
    w = grad;                   
    w(FFT_NNLS_result_zeropad == 0 & w > 0) = 0;      
    
    % Obtain the auxiliary vector g
    g = fftshift(ifft2(fft2(w).*fft_PSF));    
    
    % Estimate an optimal step factor
    step = dot(g(:),r(:))/dot(g(:),g(:));    
    
    % Update the beamforming map
    FFT_NNLS_result_zeropad = max(0,FFT_NNLS_result_zeropad - step*w);
    
end

% Remove zero padding
FFT_NNLS_result = remove_zeropad(FFT_NNLS_result_zeropad);

end

function Map_zeropad = zeropad(Map)
% To avoid wraparound effects, zero padding is used

[M,N] = size(Map);
Map_zeropad = zeros(M+N);
Map_zeropad(int64(M/2)+1:int64(M/2 + M),int64(M/2)+1:int64(M/2 + M)) = Map;

end

function Map = remove_zeropad(Map_zeropad)
% remove zeropadding

M = size(Map_zeropad,1)/2; 
Map = Map_zeropad(int64(M/2)+1:int64(M/2+M),int64(M/2)+1:int64(M/2+M));

end

function [r,g] = nnls(b,x,fft_PSF)
% Nonnegative least squares objective function
%
% Input:
%   b : beamforming map 
%   x: 	beamforming map after deconvolution
%   fft_PSF:  point spread function (PSF) after Fourier transform
%
% Output: 
%   r:  residual  
%   g:  gradient 
%

% Calculate residual 
r = fftshift(ifft2(fft2(x).*fft_PSF)) - b;

% Calculate gradient
g = fftshift(ifft2(fft2(r).*fft_PSF));

end