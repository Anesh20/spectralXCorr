% example 
clearvars
global fidm sw sfrq1H H1offset

% sample data
load sampleData;

% calling SC algorithm with
%  - spectrum range of 1.8 to 3.6ppm, 
%  - reference spectrum is the 1st transient
%  - no apodization to determine freq/phase offsets
[fidCor,outVal] = spectXcorr(fidm,[1.8 3.6],'f',0,1);



