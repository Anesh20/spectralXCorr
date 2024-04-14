function [fidCor,valOut] = spectXcorr(fid,chemicalRange, filterFlag, plotFlag)
%
%function [fidCor,valOut] = spectXcorr(fid,chemicalRange, filterFlag, plotFlag)
% Simultaneous phase and frequency estimation using cross-correlation
% in the frequency domain.  Method is termed spectral cross-correlation or SC
%
% Dinesh Deelchand, CMRR, University of Minnesota
% 24 Oct 2023
%
% INPUTS:
%     fid           - signal in time domain with dimension [np, ave]
%     chemicalRange - chemical shift range to align spectra, default is [1.8 3.5]
%     filterFlag    - apply apodization (LB=5 and GF=0.12) to  data before SC
%     plotFlag      - plot spectra and offsets (default 0)
%
% OUTPUTS:
%     fidCor - frequency and phase corrected FID
%     valOut - vector contains estimated frequency (Hz) and phase offsets (deg)
%

global sw sfrq1H H1offset

if (nargin<1)
    error('Missing input FID signal. Aborting!');
elseif (nargin==1)
    chemicalRange = [1.8  3.5];
	filterFlag = 0;
    plotFlag = 0;
elseif nargin==2
    filterFlag = 0;
    plotFlag = 0;
elseif nargin==3
    plotFlag = 0;
end

% apply LB and ZF
dw = 1/sw; t = (0:dw:dw*(length(fid)-1))';
if (filterFlag==1)
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
spect_ref = spectfft(:,1); %use first spectrum instead of mean of spectra!

% extract selected spectral region
Refchemicalranges_ppm = chemicalRange;
fmax=(sw)/2;
f=fmax:-2*fmax/(length(fidzf)-1):-fmax;
deltaFnew = sw/(length(fidzf)-1);
scale_ppm=f/(sfrq1H)+H1offset;
findval1 = find(Refchemicalranges_ppm(1)-0.1 < scale_ppm & scale_ppm < Refchemicalranges_ppm(1)+0.1, 1, 'last' );
findval2 = find(Refchemicalranges_ppm(2)-0.1 < scale_ppm & scale_ppm < Refchemicalranges_ppm(2)+0.1, 1 );
region=findval2:findval1;

ShiftCalc = zeros(1,nt);
phaseCalc = zeros(1,nt);
maxLag=round(np*(1+sifactor)/2);
tic;
for ix=1:nt
    clear c
    spect_use = (spectfft(:,ix));
    c = xcorr(spect_ref(region), spect_use(region),maxLag);
    
    %% freq shift calc
    [~,indx] = max(abs(c));
    ShiftCalcHz = (indx - (maxLag+1)) * deltaFnew;
    
    if ix==1 %relative to first scan
        ShiftRef = ShiftCalcHz;
        ShiftCalc(ix) = 0;
        indxRef=indx;
    else
        ShiftCalc(ix) = ShiftCalcHz - ShiftRef;
    end
    
    
    %% phase shift Calc
    ptsUse=5;
    if ix==1
        phzRefx = angle(c);
        %middle of max index 
        phzRef = (phzRefx(indxRef-ptsUse:indxRef+ptsUse));
        phaseCalc(ix) = 0;
    else
        phzCurx = angle(c);
        %middle of max index 
        phzCur = (phzCurx(indx-ptsUse:indx+ptsUse));
        
        %calc phase diff
        phzCal = mean(phzCur - phzRef);
        phaseCalc(ix) = rad2deg(phzCal);
    end
    
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
    curAxis=axis; axis([0 5 curAxis(3) curAxis(4)]); useAxis=axis;
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

