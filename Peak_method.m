%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the peak finding                                     %      
% interpolation method with help from matlab interp1(x,y,type)            %
% type = spline. For more info doc interp1                                % 
% type = PCHIB. Same as cubic interpolation                               %
%                                                                         %
% Input:                                                                  %
%        result(freq(ind),H1(ind),H2(ind),COHERENCE)                      %
%        omega -> xpos to peak ypos                                       %
%                                                                         %
% Output:                                                                 %
%        dratio -> damping ratios for specified data set                  %
%                                                                         %
%                                                                         %
% Made by:                                                                %
%         A group with big balls                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dratio] = Peak_method(result,omega,dratio,p)


%Find index
[index,~] = find(round(result(:,1),4) == round(omega,4));

%Set values
freq = result(:,1);
H1   = result(:,2);
H2   = result(:,3);

%% Set number of data points can't do with n = 1. Spline needs a total of 4;
n = 10;
disc = 3000; %3000 good!
type = 'PCHIB';

%Set figure
close all; figure(1); hold on
for i = 1:1:6
    
    %Initiate plotting values
    Hplot = sprintf('$H_{%d%d}/sqrt(2)$',i,p);
    OMPA  = sprintf('{a%d%d}',i,p);
    OMPA  = strcat('$','\omega_',OMPA,'$',' = ');
    OMPB  = sprintf('{b%d%d}',i,p);
    OMPB  = strcat('$','\omega_',OMPB,'$',' = ');
    
    %Set x range and y range
    xrange = [freq((index(i)-n):1:index(i):(index(i)+n))];
    yrange = [abs(H2((index(i)-n):1:index(i):(index(i)+n)))];
    
    %Create interpolated xi and yi values
    xi = linspace(freq(index(i)-n),freq(index(i)+n),disc);
    yi = interp1(xrange,yrange,xi,type);
    
    %Ploting for temporary vizualisation. Toogle of if not interested
    f1 = subplot(2,3,i); hold on
    %Ax = axes(f1)
    h3 = plot([min(xi),max(xi)],[abs(H2(index(i)))/sqrt(2),abs(H2(index(i)))/sqrt(2)],'k','LineWidth',1.2);
    h4 = plot([min(xi),max(xi)],[abs(H2(index(i))),abs(H2(index(i)))],'k','LineWidth',1.2);
    h1 = plot(xi,yi,'LineStyle','none','marker','o','markersize',5,'MarkerEdgeColor','k');
    h2 = plot(xrange,yrange,'o','markersize',12,'MarkerEdgeColor','r','LineWidth',2);   
    ylim([min(abs(yi)) 1.1*max(abs(yi))]); YA = ylim;
    xlim([min(abs(xi)) max(abs(xi))]);     XA = xlim;
    
    %legend([h2,h1,h3],{'Orignal data','PCHI data',Hplot},'FontSize',20,'Interpreter','latex')
    %annotation('doublearrow',[0.2,0.2],[abs(H2(index(i)))/sqrt(2),abs(H2(index(i)))])

    %Use formulation Inman book equation 7.13 --> |H(omega_a)|=|H(omega_b)|=|H(omega_d)|/(sqrt(2))
    H2a = abs(H2(index(i)))/sqrt(2); % 3 db downpoint
    H2b = H2a;
    
    %Determine fraction and position
    d = H2a./(yi+eps);
    ix = find(d >.985 & d <1.0125); % d >.985 & d <1.0125 Works %995 & d <1.005
    x_sol = xi(ix);
    y_sol = yi(ix);
    
    %Determine if x_sol and y_sol are on left or right side of omega(i)
    LR = x_sol-omega(i);
    SN = [];
    SP = [];
    count1 = 1;
    count2 = 1;
    for k = 1:length(LR)
        if LR(k) < 0
            SN(count1,1) = x_sol(k);
            SN(count1,2) = y_sol(k);
            SN(count1,3) = LR(k);
            SN(count1,4) = k;
            SN(count1,5) = ix(k);
            count1 = count1+1;
        elseif LR(k) > 0
            SP(count2,1) = x_sol(k);
            SP(count2,2) = y_sol(k);
            SP(count2,3) = LR(k);
            SP(count2,4) = k;
            SP(count2,5) = ix(k);
            count2 = count2+1;
        end
    end
    
    %Initialize arrays
    slope = []; xout = []; yout = [];
    %Begin gradient decent
    for k = 1:1:size(SP,1)
        %Initialize count
        count = 0;
        for j = 1:1:size(SN,1)
            %Compute gradient
            slope(j,k) = (SP(k,2)-SN(j,2))/(SP(k,1)-SN(j,1));
            %Initiate plotting
            xplot = [SN(j,1),SP(k,1)]; yplot = [SN(j,2),SP(k,2)];
            %Figure out where slope is smallest between S+ AND S-
            if abs(slope(j,k)) <= abs(slope(j-count,k))
                %Rezet and set count
                count = 0; count = count + 1;
                xout(k,:) = xplot;
                yout(k,:) = yplot;
            end
        end
    end
    
    %Creating fallping wank plots Only for vizualisation toogle of if not
    %of interest
%     for wank = 1:1:size(xout,1)
%         plot(xout(wank,:),yout(wank,:))
%     end
    
    %Initialize count
    count = 0;
    %Initialize arrays
    ed_ave = []; %outputx = []; outputy = [];
    %Finalise output
    for k = 1:1:size(yout,1)
        %Averaging
        ave = (yout(k,1)+yout(k,2))/2;
        %Euclidian distance 1D sqrt((q-p)^2) = |q-p|
        ed_ave(k) = abs(H2a-ave);
        if ed_ave(k) <= ed_ave(k-count)
            %Rezet and set count
            count = 0; count = count + 1;
            %output
            outputx(i,:) = xout(k,:);
            outputy(i,:) = yout(k,:);
            %Equation (7.14) Inmann. xi = omega_b-omega_a/(2*omega_d)
            dratio(i,:)  = (outputx(i,2) - outputx(i,1))/(2*omega(i));
        end
    end
        %Update OMPA & OMPB
        OMPA = strcat(OMPA,num2str(round(outputx(i,1),3)));
        OMPB = strcat(OMPB,num2str(round(outputx(i,2),3)));
        
        plot(outputx(i,:),outputy(i,:),'linewidt',1.25,'color','b')
        h4 = plot(outputx(i,1),outputy(i,1),'o','markersize',12,'MarkerEdgeColor','b','LineWidth',2); 
        h5 = plot(outputx(i,2),outputy(i,2),'o','markersize',12,'MarkerEdgeColor','b','LineWidth',2);
        legend([h2,h1,h4,h5],{'Orignal data','PCHI data',OMPA,OMPB},'FontSize',15,'Interpreter','latex')
end
%disp(dratio')
end