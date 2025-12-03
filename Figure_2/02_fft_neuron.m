function [peaks,majorpower]=fft_neuron(imageprofile,filename,timeinterval)

%% to run:
imageprofile=double(imageprofile);
ll=length(imageprofile);

cinfo=strrep(filename,'.tif','_output/');
[~,savename,~]=fileparts(filename);

%colorofcurve = 'r';

%scrsz = get(0,'ScreenSize');
%figure('Position',[scrsz(3)*0.5 scrsz(4)*0.5 scrsz(3)*0.3 scrsz(4)*0.3],'PaperPosition',[1 12 3.6 2]);
%plot(timeinterval:timeinterval:ll*timeinterval,imageprofile,'Color',colorofcurve,'LineWidth',0.7);
%set(gca,'XTick',0:5:round(ll*timeinterval*10)/10,'XGrid','on','YLim',[min(imageprofile) max(imageprofile)]);
%xlabel('sec');
%set(gca,'XLim',[0 round(ll*timeinterval*10)/10]);
%saveas(gca,[cinfo 'raw.eps'],'epsc');

%% substract trend component, so it has less low freq fluctuation
degree = 3;
trend = polyfit([timeinterval:timeinterval:ll*timeinterval]',imageprofile,degree); 
imageprofile = imageprofile-polyval(trend,[timeinterval:timeinterval:ll*timeinterval]');

%% use hanning filter to reduce spectral leakage
t=timeinterval:timeinterval:ll*timeinterval;
hanning = 0.5 - 0.5*cos(2*pi*t/t(end));
%figure;plot (t, imageprofile)
imageprofile= imageprofile.*hanning';
%figure;plot (t, imageprofile)

%% normalize the data
profile_ordered=sort(imageprofile);
profile_min=mean(profile_ordered(1:10));
profile_max=mean(profile_ordered(ll-1:ll));
imageprofile=(imageprofile-profile_min)/(profile_max-profile_min);

%% FFT. add padding 
n=4096;
Y = fft(imageprofile,n);

%% To ensure the operands are integers, use the ceil, fix, floor, and round functions.
P=Y.*conj(Y)/n;
power = P(1:n/2 + 1); power(2:n/2) = 2*power(2:n/2);

lengthlimit=length(power); 
nyquist = 1/2;
freq = (0:ceil(n/2))/(ceil(n/2))*nyquist;  
period=1./double(freq);                
%figure;plot(freq,power,'.');

%% to plot only the relevant part
minrepeat=5;
index_lessthan100=find(period*timeinterval<100);
firstindex=index_lessthan100(1);
%plot period
%figure('Position',[scrsz(3)*0.5 scrsz(4)*0.5 scrsz(3)*0.3 scrsz(4)*0.3],'PaperPosition',[0.6 6 2 2]);
%plot(period(minrepeat:lengthlimit-2)*timeinterval,power(minrepeat:lengthlimit-2),'Color',colorofcurve,'LineWidth',0.5);
%xlabel('Period (sec/cycle)');
%axis([0,20,0,15]) %%%% --changed the Y_Maximum_Limit_2019Sept30
%set(gca,'Xtick',[0:2:20]);
%hold on;

radius = 1;

allpeaks = zeros(length(firstindex:lengthlimit-radius));
for iii=firstindex:lengthlimit-radius
    product=1;
    for jjj=1:radius*2+1
        product=product*(power(iii)>=power((iii-radius-1+jjj)));
    end
    allpeaks(iii-firstindex+1)=product;
end


allpeaksindex=find (allpeaks==1)+firstindex-1;
areaP = trapz(freq,power);
%% find large enough peaks
hh= power(allpeaksindex)>0.02*max(power(firstindex:lengthlimit)); 
largepeaksindex=allpeaksindex(hh);
[~,IX] = sort(power(largepeaksindex),'descend');
index=largepeaksindex(IX);

%% for each large peaks, find valley so that area under the peak can be calculated
dP = gradient(power);

for j=1:min(length(index),5) % find the top 5(used to be 3) peaks that satisfy the criteria above
    majorperiod(j)=period(index(j))*timeinterval; 
    majorpower(j)=power(index(j));
    mainPeriodStr=num2str(round(period(index(j))*timeinterval*10)/10);

    % find valley/boundary for each peak
    [~,sortedValleys] = sort(abs(freq(index)-freq(index(j))));
    f1 = min(index(sortedValleys(1:2)));
    f2 = max(index(sortedValleys(1:2)));
    
    peaks(j,1:3) = [freq(index(j)) period(index(j))*timeinterval power(index(j))'];

        peaks(j,4) = 100*trapz(freq(f1:f2),power(f1:f2))/areaP;
        
    %if j==1
        %plot(period(index(j))*timeinterval,power(index(j)),'k.', 'MarkerSize',20);
        %text(period(index(j))*timeinterval+1,power(index(j)),mainPeriodStr,'FontSize', 10);
    %else
       %plot(period(index(j))*timeinterval,power(index(j)),'.','Color','k', 'MarkerSize',10);
        %text((period(index(j))-0.8)*timeinterval,power(index(j)),mainPeriodStr,'FontSize', 10,'Color','k');
    %end
    %set(gcf, 'Color','w');
    %hold on;
end

%hold off;
%saveas(gca,[cinfo 'fft_mp.eps'],'epsc');

%% to calculate relative power
peaks = sortrows(peaks,4);  
%figure;
%plot(period(minrepeat:lengthlimit-2)*timeinterval,power(minrepeat:lengthlimit-2),'+-','Color',colorofcurve,'LineWidth',0.5);  % change according to how many repeating unit are available in the image
%[zr ~]=size(peaks);
%for zzzz=1:zr
%display(['Freq=',num2str(peaks(zzzz,1)),'; Period=',num2str(peaks(zzzz,2)),'; Power=',num2str(peaks(zzzz,3)),' Percentage=',num2str(peaks(zzzz,4)),'%.']);
%end
%save([cinfo [savename '_fft.mat']],'peaks');
% close all
end