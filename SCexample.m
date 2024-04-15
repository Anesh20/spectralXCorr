% example 

clearvars


global fidm sw sfrq1H H1offset
load sampleData;

[fidCor,outVal] = spectXcorr(fidm,[1.8 3.6],0,1);