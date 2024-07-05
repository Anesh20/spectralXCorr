function [fidCor,valOut] = spectXcorr(fid,chemicalRange, ref, filterFlag, plotFlag)
%
%function [fidCor,valOut] = spectXcorr(fid,chemicalRange, ref, filterFlag, plotFlag)
% Simultaneous phase and frequency estimation using cross-correlation
% in the frequency domain.  Method is termed spectral cross-correlation or SC
%
% Dinesh Deelchand, CMRR, University of Minnesota
% 24 Oct 2023, 
% Updated 03 July 2024
%
% INPUTS:
%     fid           - signal in time domain with dimension [np, ave]
%     chemicalRange - chemical shift range to align spectra, default is [1.8 3.6]
%     refSpec       - use 1st transient (='f') or mean all spectra ('m') as reference, default is 1
%     filterFlag    - apply apodization (LB=5 and GF=0.12) to data before SC, default is off
%     plotFlag      - plot spectra and offsets, default is 0
%
% OUTPUTS:
%     fidCor - frequency and phase corrected FID
%     valOut - vector contains estimated frequency (Hz) and phase offsets (deg)
%

global sw sfrq1H H1offset

% check for input parameters 
if nargin < 1
    error('Missing input FID signal. Aborting!');
end

if nargin < 2
    chemicalRange = [1.8 3.6];
end
if nargin < 3
    ref = 'f';
end
if nargin < 4
    filterFlag = 0;
end
if nargin < 5
    plotFlag = 0;
end


% apply LB and ZF
dw = 1/sw; t = (0:dw:dw*(length(fid)-1))';
if (filterFlag==1)
    disp(' *** FID apodization applied! ***')
	LB = 5;
	GF = 0.15;	
else
	LB = 0;
	GF = 1000;
end
sifactor = 10;
[np,nt,nbCoils]=size(fid);
fidzf = complex(zeros(np*(sifactor+1),nt,nbCoils));
fidCor = complex(zeros(np,nt,nbCoils));
for ical=1:nbCoils
    for jcal=1:nt
        fidzf(:,jcal,ical) = [fid(:,jcal,ical).*exp(-t*pi*LB-t.^2/(GF^2)); zeros(np*sifactor,1)];
    end
end

% FFT and mean
spectfft = fftshift(fft(fidzf,[],1),1);
if ~strcmp(ref, 'f') && ~strcmp(ref, 'm')
    disp(' reference is out of range. Using 1st transient as reference')
    ref = 'f';
end

% reference spectrum
if (strcmp(ref,'f'))
    SRef = spectfft(:,1); %use first spectrum 
    disp(' Using 1st transient as reference')
elseif (strcmp(ref,'m'))
    SRef = mean(spectfft,2)/size(spectfft,2); %use mean spectra 
    disp(' Using mean spectrum as reference')
end

% extract selected spectral region
Refchemicalranges_ppm = chemicalRange;
fmax=(sw)/2;
f=fmax:-2*fmax/(length(fidzf)-1):-fmax;
deltaFnew = sw/(length(fidzf)-1);
scale_ppm=f/(sfrq1H)+H1offset;
findval1 = find(Refchemicalranges_ppm(1)-0.1 < scale_ppm & scale_ppm < Refchemicalranges_ppm(1)+0.1, 1, 'last' );
findval2 = find(Refchemicalranges_ppm(2)-0.1 < scale_ppm & scale_ppm < Refchemicalranges_ppm(2)+0.1, 1 );
region = findval2:findval1;

ShiftCalc = zeros(1,nt);
phaseCalc = zeros(1,nt);
maxLag=round(np*(1+sifactor)/2);

tic;

%% reference scan correlation 
CRef = xcorr(SRef(region), SRef(region),maxLag);
[~,indx] = max(abs(CRef)); %freq ref
ShiftRef = (indx - (maxLag+1)) * deltaFnew;
indxRef = indx;
phzRefx = angle(CRef); %phase ref
ptsUse=5; %on each side
phzRef = (phzRefx(indxRef-ptsUse:indxRef+ptsUse));
 
for ix=1:nt
    clear Sn Cn
    Sn = (spectfft(:,ix));
    Cn = xcorr(SRef(region), Sn(region),maxLag);
    
    %% freq shift calc
    [~,indx] = max(abs(Cn));
    ShiftCalcHz = (indx - (maxLag+1)) * deltaFnew;
    ShiftCalc(ix) = ShiftCalcHz - ShiftRef;
    
    %% phase shift Calc
    phzCurx = angle(Cn);
    %middle of max index
    phzCur = (phzCurx(indx-ptsUse:indx+ptsUse));
    
    %calc phase diff
    phzCal = mean(phzCur - phzRef);
    phaseCalc(ix) = rad2deg(phzCal);
    
	%% Freq and phase corrected FID
    fidCor(:,ix) = fid(:,ix).*exp(1i*2*pi*ShiftCalc(ix).*t).*exp(1i*deg2rad(phaseCalc(ix)));
end
elapsed_time = toc * 1000;
fprintf('SC Time taken: %.2f ms\n', elapsed_time);

% output
valOut = [ShiftCalc; phaseCalc]';

if plotFlag
    f=fmax:-2*fmax/(length(fid)-1):-fmax;
    scale_ppmOrig = f/(sfrq1H)+H1offset;
    figure, clf
    
    spectfftOrig = fftshift(fft(fid,[],1),1);
    subplot(221), plot(scale_ppmOrig,real(spectfftOrig)); title('Original data');
    set(gca,'xdir','reverse')
    curAxis=axis; axis([0.5 4.5 curAxis(3) curAxis(4)]); useAxis=axis;
    xlabel('Chemical shift (ppm)')
    %average spectrum
    hold on, plot(scale_ppmOrig,mean(real(spectfftOrig),2),'k','linewidth',2); 
    
    spectfftCor = fftshift(fft(fidCor,[],1),1);
    subplot(222), plot(scale_ppmOrig,real(spectfftCor)); title('Corrected data'); set(gca,'xdir','reverse')
    axis(useAxis);
    xlabel('Chemical shift (ppm)')
    hold on, plot(scale_ppmOrig,mean(real(spectfftCor),2),'k','linewidth',2); 
    
    subplot(223), plot(-ShiftCalc); title('Frequency shift (Hz)');
    xlabel('Scan number')
    
    subplot(224), plot(-phaseCalc), title('Phase offset (deg)');
    xlabel('Scan number')
end

return

