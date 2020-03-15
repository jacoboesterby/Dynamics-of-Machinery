

%This file contains data found from visual inspection of the 6x2 different
%datafiles


%% Find frequency peak for 6 different data files
omega1 = [0.84,1.22,1.77,3.64,5.44,7.26];   %Frequencies found from visual inspection of dataset 1.
omega2 = [0.85,1.22,1.76,3.61,5.52,7.20];   %Frequencies found from visual inspection of dataset 2.
omega3 = [0.84,1.22,1.76,3.63,5.51,7.18];   %Frequencies found from visual inspection of dataset 3.
omega4 = [0.84,1.22,1.76,3.61,5.53,7.24];   %Frequencies found from visual inspection of dataset 4.
omega5 = [0.82,1.18,1.76,3.63,5.52,7.25];   %Frequencies found from visual inspection of dataset 5.
omega6 = [0.84,1.21,1.73,3.63,5.49,7.25];   %Frequencies found from visual inspection of dataset 6.
omega_nat = [omega1;omega2;omega3;omega4;omega5;omega6];
%% Find transfer functions
H11 = [0.1205,0.1574,0.689,4.561,56.69,7.062];  %Transfer functions found from visual inspection of dataset 1.
H21 = [0.5328,0.4706,6.692,56.15,7.827,35.3];   %Transfer functions found from visual inspection of dataset 1.
H31 = [2.122,1.029,5.246,27.01,29.73,27.9];     %Transfer functions found from visual inspection of dataset 1.
H41 = [1.777,2.689,3.959,28.74,13.42,6.163];    %Transfer functions found from visual inspection of dataset 1.        
H51 = [11.88,7.255,6.647,4.474,1.48,0.2674];    %Transfer functions found from visual inspection of dataset 1.
H61 = [3.922,2.588,10.24,4.284,0.8499,0.3213];  %Transfer functions found from visual inspection of dataset 1.
H = [H11;H21;H31;H41;H51;H61];