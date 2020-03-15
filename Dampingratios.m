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
%       Modified by a group with big balls              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Run data
run data.m

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


%% IMPORTANT PARAMETERS OBTAINED FROM THE ACQUISITION SYSTEM
N        = 5000; % number of points in each block of the total signal
Fs       = 50;   % sampling frequency used while acquiring the signals [Hz]
WIN = N;	     % number of the points used by the HANNING (default) window
%WIN(1:N) = 1;	 % number of the points used by the HANNING (default)window
OVLP     = 0;	 % "overlap" factor among the blocks (0 until 1)
SENS_X   = 1/1;  %/0.316;   % sensibility of the input signal [N/V]
SENS_Y   = 1/1;  % sensibility of the output signal [m/V]


%% LOADING THE INPUT AND OUTPUT SIGNALS FROM FILES ~ (Modification loop)
for j = 1:1:6
    force = strcat('force',num2str(j),'.txt');
    acc   = strcat('acceleration',num2str(j),'.txt'); 
    
    data(j).x = load(force)*SENS_X;
    data(j).y = load(acc)*SENS_Y;
    
    data(j).x_org = data(j).x;
    data(j).y_org = data(j).y;
end 

%% FILTER DESIGN 
  [Nz,Dz]=butter(5,0.99);   % filter (Butterworth) 
                            % 5th order filter with 
                            % a cut-off frequency of 0.99*(Fs/2)= 24.75 Hz  
                           
%% FILTERING OF INPUT AND OUTPUT SIGNALS  
for j = 1:1:6  
   data(j).x = filter(Nz,Dz,data(j).x); % input signal after filtering
   data(j).y = filter(Nz,Dz,data(j).y); % output signal after filtering  
end

for j = 1:1:6    
   %Power spectral densities
   data(j).PXY = cpsd(data(j).x,data(j).y,WIN,OVLP,N,Fs);
   data(j).PYX = cpsd(data(j).y,data(j).x,WIN,OVLP,N,Fs); 
   
   %Cross spectral densities
   data(j).PXX = pwelch(data(j).x,WIN,OVLP,N,Fs);
   data(j).PYY = pwelch(data(j).y,WIN,OVLP,N,Fs);
   

    % FRF ESTIMATION USING THE ESTIMATORS H1 AND H2 
    data(j).H1y = (data(j).PXY)./(data(j).PXX);
    data(j).H2y = (data(j).PYY)./(data(j).PYX);
    
    % Frequency & Coherence
    [data(j).Cy,data(j).F] = mscohere(data(j).x,data(j).y,WIN,OVLP,N,Fs);
    
    % GRHPHICS OF THE PROCESSED SIGNALS IN TIME AND FREQUENCY DOMAINS
    data(j).t = [1:size(data(j).x)]/Fs; %time
    idx = [1:N/6];                      % N/2 -> 25 Hz  or  N/4 -> 12.5 Hz  
    
    % Calculating the phase based on the function ACOS.
    for cont1=1:max(size(idx))
        data(j).fase1y(cont1) = acos( real(data(j).H1y(cont1)) / sqrt( (real(data(j).H1y(cont1)))^2+(imag(data(j).H1y(cont1)))^2) );
        data(j).fase2y(cont1) = acos( real(data(j).H2y(cont1)) / sqrt( (real(data(j).H2y(cont1)))^2+(imag(data(j).H2y(cont1)))^2) );
    end   
end

%_____________________________________________________
% Printing the final results in the form of table: 
% FREQ. [Hz]    H1 [(m/s^2)/N]     H2 [(m/s^2)/N]    COHERENCE
result1 = [data(1).F(idx), data(1).H1y(idx), data(1).H2y(idx), data(1).Cy(idx)];
result2 = [data(2).F(idx), data(2).H1y(idx), data(2).H2y(idx), data(2).Cy(idx)];
result3 = [data(3).F(idx), data(3).H1y(idx), data(3).H2y(idx), data(3).Cy(idx)];
result4 = [data(4).F(idx), data(4).H1y(idx), data(4).H2y(idx), data(4).Cy(idx)];
result5 = [data(5).F(idx), data(5).H1y(idx), data(5).H2y(idx), data(5).Cy(idx)];
result6 = [data(6).F(idx), data(6).H1y(idx), data(6).H2y(idx), data(6).Cy(idx)];

%% Plotting 
% Graphics 1 - Analysis in the time domain
% Graphics 2 - Analysis in the frequency domain
% Graphics 3 - Analysis in the frequency domain - H2
close all
PP = 2;

exp_plotfunc(data(PP).x,      data(PP).y,...
             data(PP).t,      data(PP).F,...
             data(PP).H1y,    data(PP).H2y,...
             data(PP).fase1y, data(PP).fase2y,...
             data(PP).Cy,     idx)

%% Peak method ~ use omega(p) where p denotes accelerated mass
for i = 1:1:6
    p = i; xi1 = []; result = strcat('result',num2str(i));   
    xi1 = Peak_method(eval(result),omega_nat(i,:),xi1,p);
    %Put results into array
    xiPM(i,:) =xi1;
end
         
%% LS method ~ use omega(p) where p denoted accelerated mass
% n is number of points on either side of natfreq, Ilmar recomeds 3.
for i = 1:1:6
    p = i; xi2 = []; n = 4; result = strcat('result',num2str(i));        
    xi2 = LS_method(eval(result),omega_nat(1,:),xi2,p,n);
    %Put results into array
    xiLS(i,:) = xi2;
end

%Averaging the damping ratios for each mass
xiPMA = [mean(xiPM(:,1)),mean(xiPM(:,2)),mean(xiPM(:,3)),...
         mean(xiPM(:,4)),mean(xiPM(:,5)),mean(xiPM(:,6))];
xiLSA = [mean(xiLS(:,1)),mean(xiLS(:,2)),mean(xiLS(:,3)),...
         mean(xiLS(:,4)),mean(xiLS(:,5)),mean(xiLS(:,6))];
    
%Median of the damping for each mass
xiPMM = [median(xiPM(:,1)),median(xiPM(:,2)),median(xiPM(:,3)),...
         median(xiPM(:,4)),median(xiPM(:,5)),median(xiPM(:,6))];    
xiLSM = [median(xiLS(:,1)),median(xiLS(:,2)),median(xiLS(:,3)),...
         median(xiLS(:,4)),median(xiLS(:,5)),median(xiLS(:,6))];    

%Display damping ratios
disp('Damping ratios Peak method ~ each row reprenst sample data on mass n')
disp(xiPM)
disp('Damping ratios Peak method ~ Averaged and median')
disp(xiPMA) 
disp(xiPMM)
disp('Damping ratios LS method ~ each row reprenst sample data on mass n')
disp(xiLS)
disp('Damping ratios LS method ~ Averaged and median')
disp(xiLSA) 
disp(xiLSM)