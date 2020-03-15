%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the Least square method                              %
%                                                                         %
% Input:                                                                  %
%        result(freq(ind),H1(ind),H2(ind),COHERENCE)                      %
%        omega -> xpos to peak ypos                                       %
%                                                                         %
% Output:                                                                 %
%        dratio  -> damping ratios for specified data set                 %
%        natfreq -> calculated natural frequencies (not toogled on)       %
%                                                                         %
%                                                                         %
% Made by:                                                                %
%         A group with big balls                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [dratio] = LS_method(result,omega,dratio,p,n)

% %Normalize n
% nind = floor(n/2);

%Find index
[index,~] = find(round(result(:,1),4) == round(omega,4));

%Set values
freq = result(:,1);
H1   = result(:,2);
H2   = result(:,3);

%Bulding matrix vector system Ax = b
for i = 1:1:6
    
    
    %Find associated values
    xrange = [freq((index(i)-n):1:index(i):(index(i)+n))];
    yrange = [H2((index(i)-n):1:index(i):(index(i)+n))];
    
    %Creating size
    sx = length(xrange);
    
    %Prealocate arrays
    A  = ones(sx,2);  b  = zeros(sx,1);
    At = zeros(sx,1); bt = zeros(sx,1);
    
    %For obtaining {m,k}^T
    A(:,1)  = -1*xrange.^2;
    b(:,1)  = real(xrange.^2./(yrange));
    mk      = A\b;     %(transpose(A)*A)^(-1)*transpose(A)*b;
    
    %For obtaining {d}
    At(:,1) = xrange;
    bt(:,1) = imag(xrange.^2./(yrange));
    d       = At\bt;   %(transpose(At)*At)^(-1)*transpose(At)*bt;
    
    %Calculating natural frequencies and damping ratios
    dratio(i)  = abs(d/(2*sqrt(mk(2)*mk(1)))); %Calculating damping ratio
    natfreq(i) = sqrt(mk(2)/mk(1));            %Calculating natural frequency
end

% %Printing results
% vis1 = sprintf('Natural frequencies for data %d',p);
% vis2 = sprintf('Damping ratios for data %d and number of points %d',p,n);
% disp(vis1)
% disp(natfreq)
% disp(vis2)
% disp(dratio)

end
