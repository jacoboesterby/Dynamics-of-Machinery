function AccTimeDomain(accel,force,lab)

%
% IMPORTANT: Remember to change the extention of this file from *.txt to *.m 
%            before running Matlab! 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 41614 DYNAMICS OF MACHINES                            %
% MEK - DEPARTMENT OF MECHANICAL ENGINEERING            %
% DTU - TECHNICAL UNIVERSITY OF DENMARK                 %
%                                                       %
%               Copenhagen, March 7th,    2008          %
%                                                       %
%                                                       %
% PROJECT 1 - HELP-ROTINE TO PLOT THE ACCELERATION      %
%             SIGNAL IN TIME DOMAIN                     %
%                                                       %
%                            Ilmar Ferreira Santos      %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% IMPORTANT PARAMETERS OBTAINED FROM THE ACQUISITION SYSTEM

N      = 6000; % number of points in each block of the total signal
Fs     = 50;   % sampling frequency used while acquiring the signals [Hz]
SENS_Y = 1;    % sensibility of the output signal [m/s^2/V]
SENS_F = 1;    % sensibility of the input signal [m/s^2/N]

%_______________________________________________________________________  
% LOADING THE ACCELERATION SIGNALS FROM FILES 
% load acceleration.txt;    % loading the file acceleration.txt wiht the output signal 
% y = acceleration*SENS_Y;  % output signal (before filtering) multiplied 
                            % by the sensor sensibility (acceleration sensor) 
  
                            
% CASE 1 : TRANSIENT ANALYSIS - MODE SHAPE 1
                               
y = accel*SENS_Y; 
                               
f = force*SENS_F;   

% % CASE 2 : TRANSIENT ANALYSIS - MODE SHAPE 2
% load acceleration_transient_mode2.txt;                                
% y = acceleration_transient_mode2*SENS_Y; 
% load force_transient_mode2.txt;                                
% f = force_transient_mode2*SENS_F; 
% 
% % CASE 3 : TRANSIENT ANALYSIS - MODE SHAPE 3
% load acceleration_transient_mode3.txt;                                
% y = acceleration_transient_mode3*SENS_Y; 
% load force_transient_mode3.txt;                                
% f = force_transient_mode3*SENS_F; 
%  
% % CASE 4 : TRANSIENT ANALYSIS - MODE SHAPE 4
% load acceleration_transient_mode4.txt;                                
% y = acceleration_transient_mode4*SENS_Y; 
% load force_transient_mode4.txt;                                
% f = force_transient_mode4*SENS_F; 
%  
% % CASE 5 : TRANSIENT ANALYSIS - MODE SHAPE 5
% load acceleration_transient_mode5.txt;                                
% y = acceleration_transient_mode5*SENS_Y; 
% load force_transient_mode5.txt;                                
% f = force_transient_mode5*SENS_F; 
% 
% % CASE 6 : TRANSIENT ANALYSIS - MODE SHAPE 6
% load acceleration_transient_mode6.txt;                                
% y = acceleration_transient_mode6*SENS_Y; 
% load force_transient_mode6.txt;                                
% f = force_transient_mode6*SENS_F; 
                   
%_______________________________________________________________________                            

% GRHPHICS OF THE SIGNAL IN TIME DOMAIN
  
t=[1:size(y)]/Fs;               % time vector 
Freq=[0:Fs/N:Fs/2];              % frequency vector
y_fft=abs(fft( y(1:N) ) )/N;    % fft of the acceleration signal
f_fft=abs(fft(f(1:N)))/N;       % fft of the force signal

figure
plot(t(1:N),y(1:N),'r-','LineWidth',1.5); 
grid
set(gca,'FontAngle','oblique','FontSize',14)
title(strcat('Acceleration of Mass - ',lab))
ylabel('acc [m/s^2]')
xlabel('time [s]')

figure
plot(t(1:N),f(1:N),'r-','LineWidth',1.5); 
grid
set(gca,'FontAngle','oblique','FontSize',14)
title(strcat('Force on Mass 1 - ',lab))
ylabel('force [N]')
xlabel('time [s]')

%_______________________________________________________________________                            

% GRHPHICS OF THE SIGNAL IN FREQUENCY DOMAIN

figure
subplot(2,1,1); plot(t(1:N),f(1:N),'r-','LineWidth',1.5); 
grid
set(gca,'FontAngle','oblique','FontSize',14)
title(strcat('(a) Force in Time Domain - ',lab))
ylabel('force [N]')
xlabel('time [s]')
subplot(2,1,2); plot(Freq(1:N/5),f_fft(1:N/5),'r-','LineWidth',1.5);
grid
set(gca,'FontAngle','oblique','FontSize',14)
title(strcat('(b) Force in Frequency Domain - ',lab))
ylabel('FFT(force) [N]')
xlabel('Freq [Hz]')

figure
subplot(2,1,1); plot(t(1:N),y(1:N),'r-','LineWidth',1.5); 
grid
set(gca,'FontAngle','oblique','FontSize',14)
title(strcat('(c) Acceleration in Time Domain - ',lab))
ylabel('acc [m/s^2]')
xlabel('time [s]')
subplot(2,1,2); plot(Freq(1:N/5),y_fft(1:N/5),'r-','LineWidth',1.5);
grid
set(gca,'FontAngle','oblique','FontSize',14)
title(strcat('(d) Acceleration in Frequency Domain - ',lab))
ylabel('FFT(acc) [m/s^2]')
xlabel('Freq [Hz]')



