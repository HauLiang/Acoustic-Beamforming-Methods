function [DAS_result, PSF, hn, CSM] = DAS(N,z0,f,phi,rn,source,SNR)
%
% This code implements the DAS algorithm
%
% More information about DAS can be found in the paper:
%    Van Veen, Barry D and Buckley, Kevin M, 
%    "Beamforming: A versatile approach to spatial filtering", 
%    IEEE assp magazine, 1988.
%
%
% Inputs:
%    N:  number of grid points in each dim 
%    z0: source distance
%    f:  sampling frequency
%    rn: coordinates of the microphone array
%    source: x,y position of sources
%    SNR:  signal-to-noise ratio
%
% Outputs:
%    DAS_result:  beamforming map, obtained by DAS
%    PSF:  point spread function (PSF)
%    hn:   steering vector
%    CSM: cross-spectrum matrix (CSM)
%
% Author: Hao Liang 
% Last modified by: 21/09/07
%

% Number of microphones in array
N_mic = size(rn,1);

% Parameters
c = 343;         % Speed of sound
omega = 2*pi*f;  % Angular frequency
Np = round(N*1.3);     % PSF grid size

% Parameters initialization 
dj = zeros(Np,Np,N_mic); PSF = zeros(Np,Np);
hn = zeros(Np,Np,N_mic); gn = zeros(Np,Np,N_mic); 
DAS_result = zeros(N,N);

% Scan plane
L = 2*z0*tand(phi);            
x = [-L/2 L/2];    
scan_range = linspace(x(1),x(2),N);
[X,Y] = meshgrid(scan_range);

% PSF grid
rx_psf = linspace(x(1),x(2),Np);
[Xp,Yp] = meshgrid(rx_psf);

% |d0|: Distance from (0,0,0) to all mesh points
d0 = sqrt(Xp.^2 + Yp.^2 + z0^2);    

% Distance dj from each microphone to each grid point
for n = 1:N_mic
    dj(:,:,n) = sqrt((Xp-rn(n,1)).^2+(Yp-rn(n,2)).^2 + z0^2);
    hn(:,:,n) = (dj(:,:,n)./d0).*exp(1j*omega.*dj(:,:,n)./c);
    gn(:,:,n) = (d0./dj(:,:,n)).*exp(1j*omega.*dj(:,:,n)./c);
end

% Point spread function (PSF)
% Constructed by the DAS beamformer response of a unit strength point
% source in the middel of the grid
ind = round(Np/2);
e_unit = squeeze(gn(ind,ind,:));

for ii = 1:length(Xp)
    for jj = 1:length(Yp)
        PSF(ii,jj) = dot(squeeze(hn(ii,jj,:)),e_unit);
    end
end

PSF = rot90(abs(PSF).^2/N_mic^2);      % Normalize PSF
if N ~= Np
    PSF = interp2(Xp,Yp,PSF,X,Y);
end

% Cross spectral matrix (CSM)
dj = zeros(N,N,N_mic);
hn = zeros(N,N,N_mic);
gn = zeros(N,N,N_mic);

% Distance and steering vectors
d0 = sqrt(X.^2 + Y.^2 + z0^2);    % |d0|: Distance from (0,0,0) to all mesh points
for n = 1:N_mic
    dj(:,:,n) = sqrt((X-rn(n,1)).^2+(Y-rn(n,2)).^2 + z0^2);
    hn(:,:,n) = (dj(:,:,n)./d0).*exp(1j*omega.*dj(:,:,n)./c);
    gn(:,:,n) = (d0./dj(:,:,n)).*exp(1j*omega.*dj(:,:,n)./c)+10^(-SNR/10)*(rand(N,N)+1j*rand(N,N));
end

% Create CSM
q = zeros(N,N);
CSM = zeros(N_mic,N_mic);

for k = 1:size(source,1)
    q(source(k,2),source(k,1)) = 1;
    CSM = CSM + squeeze(gn(source(k,2),source(k,1),:))*squeeze(gn(source(k,2),source(k,1),:))';
end

% Diagonal removal technique
CSM(logical(eye(size(CSM)))) = 0;

% DAS imaging
for ii = 1:length(X)
    for jj = 1:length(Y)
        e = squeeze(hn(ii,jj,:));
        DAS_result(ii,jj) = dot(e,CSM*e)/(N_mic^2-N_mic);
    end
end

end
