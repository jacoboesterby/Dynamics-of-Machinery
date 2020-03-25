function GenerateFRF(accel,force,lab)

%
% IMPORTANT: rememder to rename the file frf_general_rio_modified.txt to frf_general_rio_modified.m before running! 
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 41514 DYNAMICS OF MACHINARY                           %
% MEK - DEPARTMENT OF MECHANICAL ENGINEERING            %
% DTU - TECHNICAL UNIVERSITY OF DENMARK                 %
%                                                       %
%              Copenhagen, February 1st,    2018        %
%                                                       %
%                                                       %
% PROJECT 1 - ROTINE TO ESTIMATE FRF AND COHERENCE      %
%             USING MATLAB SIGNAL ANALYSIS TOOLBOX      %
%                                                       %
%                            Ilmar Ferreira Santos      %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%close all


% OBSERVATION: The whole signal in time domain is divided 
% into blocks of size "N" (points) 
% Hanning window is used 
% Averaging is also done, and it is depending on the number 
% of points of signal in time domain, on the number of points 
% in the blocks and finally on the overlap factor.
% For faciliate the understanding of the code, the nomenclature 
% used is:  
% x for the input signal  (normally force signal)
% y for the output signal (normally acceleration)
% z for the output signal (normally acceleration)


% IMPORTANT PARAMETERS OBTAINED FROM THE ACQUISITION SYSTEM

N        = 5000; % number of points in each block of the total signal
Fs       = 50;   % sampling frequency used while acquiring the signals [Hz]
WIN = N;	     % number of the points used by the HANNING (default) window
%WIN(1:N) = 1;	 % number of the points used by the HANNING (default)window
OVLP     = 0;	 % "overlap" factor among the blocks (0 until 1)
SENS_X   = 1/1;  %/0.316;   % sensibility of the input signal [N/V]
SENS_Y   = 1/1;  % sensibility of the output signal [m/V]


% LOADING THE INPUT AND OUTPUT SIGNALS FROM FILES 
                                                
%load force1.txt;          
%load acceleration1.txt;     
x = force*SENS_X;                                 
y = accel*SENS_Y;                            
x_original = x;
y_original = y;

% FILTER DESIGN 

% BUTTER Butterworth digital and analog filter design.
% [B,A] = BUTTER(N,Wn) designs an Nth order lowpass digital
% Butterworth filter and returns the filter coefficients in length 
% N+1 vectors B (numerator) and A (denominator). The coefficients 
% are listed in descending powers of z. The cutoff frequency 
% Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to 
% half the sample rate.

  [Nz,Dz]=butter(5,0.99);   % filter (Butterworth) 
                            % 5th order filter with 
                            % a cut-off frequency of 0.99*(Fs/2)= 24.75 Hz  
                           
% FILTERING OF INPUT AND OUTPUT SIGNALS  

% FILTER One-dimensional digital filter.
% Y = FILTER(B,A,X) filters the data in vector X with the
% filter described by vectors A and B to create the filtered
% data Y.  The filter is a "Direct Form II Transposed"
% implementation of the standard difference equation:
% 
% a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                       - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% 
% If a(1) is not equal to 1, FILTER normalizes the filter
% coefficients by a(1). 
% FILTER always operates along the first non-singleton dimension,
% namely dimension 1 for column vectors and non-trivial matrices,
% and dimension 2 for row vectors.
                           
   x = filter(Nz,Dz,x);    % input signal after filtering
   y = filter(Nz,Dz,y);  	% output signal after filtering 


% POWER SPECTRAL DENSITY & CROSS SPECTRAL DENSITY 

% PSD Power Spectral Density estimate.
% Pxx = PSD(X,NFFT,Fs,WINDOW) estimates the Power Spectral Density of 
% a discrete-time signal vector X using Welch's averaged, modified
% periodogram method.     
% X is divided into overlapping sections, each of which is detrended 
% (according to the detrending flag, if specified), then windowed by 
% the WINDOW parameter, then zero-padded to length NFFT.  The magnitude 
% squared of the length NFFT DFTs of the sections are averaged to form
% Pxx.  Pxx is length NFFT/2+1 for NFFT even, (NFFT+1)/2 for NFFT odd,
% or NFFT if the signal X is complex.  If you specify a scalar for 
% WINDOW, a Hanning window of that length is used.  Fs is the sampling
% frequency which doesn't affect the spectrum estimate but is used 
% for scaling the X-axis of the plots.

% CSD Cross Spectral Density estimate.
% Pxy = CSD(X,Y,NFFT,Fs,WINDOW) estimates the Cross Spectral Density of 
% signal vectors X and Y using Welch's averaged periodogram method.  X and
% Y are divided into overlapping sections, each of which is detrended, 
% then windowed by the WINDOW parameter, then zero-padded to length NFFT.
% The products of the length NFFT DFTs of the sections of X and Y are 
% averaged to form Pxy.  Pxy is length NFFT/2+1 for NFFT even, (NFFT+1)/2
% for NFFT odd, or NFFT if the either X or Y is complex.  If you specify 
% a scalar for WINDOW, a Hanning window of that length is used.  Fs is the
% sampling frequency which doesn't effect the cross spectrum estimate 
% but is used for scaling of plots.
    
%  PXY = csd(x,y,N,Fs,WIN,OVLP);
%  PYX = csd(y,x,N,Fs,WIN,OVLP);
  PXY = cpsd(x,y,WIN,OVLP,N,Fs);
  PYX = cpsd(y,x,WIN,OVLP,N,Fs);
%  PXX = psd(x,N,Fs,WIN,OVLP);
%  PYY = psd(y,N,Fs,WIN,OVLP);
  PXX = pwelch(x,WIN,OVLP,N,Fs);
  PYY = pwelch(y,WIN,OVLP,N,Fs);

% FRF ESTIMATION USING THE ESTIMATORS H1 AND H2 

 H1y = PXY./PXX;
 H2y = PYY./PYX;

% CALCULATION OF THE COHERENCE FUNCTION

% COHERE Coherence function estimate.
% Cxy = COHERE(X,Y,NFFT,Fs,WINDOW) estimates the coherence of X and Y
% using Welch's averaged periodogram method.  Coherence is a function
% of frequency with values between 0 and 1 that indicate how well the
% input X corresponds to the output Y at each frequency.  X and Y are 
% divided into overlapping sections, each of which is detrended, then 
% windowed by the WINDOW parameter, then zero-padded to length NFFT.  
% The magnitude squared of the length NFFT DFTs of the sections of X and 
% the sections of Y are averaged to form Pxx and Pyy, the Power Spectral
% Densities of X and Y respectively. The products of the length NFFT DFTs
% of the sections of X and Y are averaged to form Pxy, the Cross Spectral
% Density of X and Y. The coherence Cxy is given by
%
%           Cxy = (abs(Pxy).^2)./(Pxx.*Pyy)
%
% Cxy has length NFFT/2+1 for NFFT even, (NFFT+1)/2 for NFFT odd, or NFFT
% if X or Y is complex. If you specify a scalar for WINDOW, a Hanning 
% window of that length is used.  Fs is the sampling frequency which does
% not effect the cross spectrum estimate but is used for scaling of plots.
% 
% [Cxy,F] = COHERE(X,Y,NFFT,Fs,WINDOW,NOVERLAP) returns a vector of freq-
% uencies the same size as Cxy at which the coherence is computed, and 
% overlaps the sections of X and Y by NOVERLAP samples.
    
%  [Cy,F] = cohere(x,y,N,Fs,WIN,OVLP);

[Cy,F] = mscohere(x,y,WIN,OVLP,N,Fs);
  
% GRHPHICS OF THE PROCESSED SIGNALS IN TIME AND FREQUENCY DOMAINS
  
%_____________________________________________________  
idx = [1:N/6];     % N/2 -> 25 Hz  or  N/4 -> 12.5 Hz   
t=[1:size(x)]/Fs;  % time vector 

% Calculating the phase based on the function ACOS.
for cont1=1:max(size(idx))
fase1y(cont1)=acos( real(H1y(cont1)) / sqrt( (real(H1y(cont1)))^2+(imag(H1y(cont1)))^2) );
fase2y(cont1)=acos( real(H2y(cont1)) / sqrt( (real(H2y(cont1)))^2+(imag(H2y(cont1)))^2) );
end

%_____________________________________________________
% Graphics 1 - Analysis in the time domain
figure(1)

subplot(2,1,1); 
plot(t,y,'-');
hold on
grid
title(strcat('(a) Acceleration in Time Domain - ',lab))
ylabel('(a) acc [m/s^2]')
xlabel('time [s]')
axis([0 max(t) -5 5])
subplot(2,1,2); 
plot(t,x,'-');
hold on
grid
title(strcat('(b) Force Signal in Time Domain - ',lab))
ylabel('(b) force [N]')
xlabel('time [s]')
axis([0 max(t) -5 5])


%_____________________________________________________
% Graphics 2 - Analysis in the frequency domain
figure(2)

subplot(3,1,1); 
plot(F(idx),Cy(idx),'-','LineWidth',1.5);
hold on
axis([0 max(F(idx)) 0 1.1])
title(strcat('Frequency domain - ',lab))
set(gca,'FontAngle','oblique','FontSize',14)
ylabel('Coherence ')
grid

subplot(3,1,2); 
plot(F(idx),abs(H1y(idx)),'b-',F(idx),abs(H2y(idx)),'r-','LineWidth',1.5);
hold on
axis([0 max(F(idx)) 0  1.1*max(abs(H2y(idx)))])
grid
set(gca,'FontAngle','oblique','FontSize',14)
ylabel('FRF [(m/s^2)/N] ')
legend('H1','H2')
subplot(3,1,3);
plot(F(idx),fase1y(idx)*180/pi,'b-',F(idx),fase2y(idx)*180/pi,'r-','LineWidth',1.5);
hold on
axis([0 max(F(idx)) 0 200])
grid
set(gca,'FontAngle','oblique','FontSize',14,'YTick',[-180:90:180],'YTickLabel',[-180:90:180])
set(gca,'YTick',[-180:90:180],'YTickLabel',[-180:90:180])
axis([0 max(F(idx)) 0 200])
xlabel('Frequency [Hz]')
ylabel('Phase [^o]')
legend('H1','H2')

%_____________________________________________________
% Graphics 3 - Analysis in the frequency domain - H2
figure(3)

subplot(2,1,1); 
semilogy(F(idx),abs(H2y(idx)),'-','LineWidth',1.5);
hold on
axis([0.7 max(F(idx)) 0  1.1*max(abs(H2y(idx)))])
grid
title(strcat('Frequency domain - ',lab))
set(gca,'FontAngle','oblique','FontSize',14)
ylabel('FRF [(m/s^2)/N] ')
subplot(2,1,2); 
plot(F(idx),fase2y(idx)*180/pi,'-','LineWidth',1.5);
hold on
grid
set(gca,'FontAngle','oblique','FontSize',14,'YTick',[-180:90:180],'YTickLabel',[-180:90:180])
set(gca,'YTick',[-180:90:180],'YTickLabel',[-180:90:180])
axis([0.7 max(F(idx)) 0 200])
xlabel('Frequency [Hz]')
ylabel('Phase [^o]')

%_____________________________________________________
% Printing the final results in the form of table: 
% FREQ. [Hz]    H1 [(m/s^2)/N]     H2 [(m/s^2)/N]    COHERENCE

 RESULT = [F(idx), H1y(idx), H2y(idx), Cy(idx)];