function [] = IntensityAnalysis(filename,ROImatrix,filetime,sigma_threshold,colors)
data=mean(ROImatrix');
[num_record,data_len]=size(data);
timeinterval = filetime/data_len;
amplitudes = cell(num_record,1);

[~,savename,~]=fileparts(filename);
cd(savename);

cd('../');
for j=1:num_record
    spike_data = data(j,:);
    sigma = std(data(j,:)); 
    threshold = mean(data(j,:)) + sigma_threshold * sigma;
    spike_index = spike_detection(spike_data, threshold);
    amplitude = data(j, spike_index);
    timestamp{j} = sparse(j,spike_index,amplitude);
    amplitudes{j} = amplitude;
end
cd(savename);

mid=round(data_len/2);
for m = 1:num_record
    min1 = min(data(m,1:100));
    min2 = min(data(m,(mid-50):(mid+50)));
    min3 = min(data(m,(data_len-100):(data_len)));
    mins = [min1 min2 min3];
    F0(m,1) = mean(mins);
    amplitudes_dff{m} = (amplitudes{m} - F0(m,1))./F0(m,1);
    mean_Sync_amplitude_dFF(m,1) = mean((amplitudes{m} - F0(m,1))./F0(m,1));
end
std_Sync_dff = std(amplitudes_dff{1,:});

[~, num_SynPeaks] = size(amplitudes{1,1});
Sync_rate= num_SynPeaks / (filetime/60);

figure;
hold on;
for i = 1:num_record
    plot((num_record - i) + timestamp{i}(i,:), '-', 'linewidth', 2,'Color','k');
end
xlim([0 data_len]);
set(gca, 'YTick', [], 'ycolor', 'k');
ylabel('Synchronous spikes','Fontname','Arial', 'FontSize', 10);
xlabel('Frame','Fontname','Arial', 'FontSize', 10);
hold off;
print('-dtiff',strcat(savename, '_synchspikes.png'));

% Network Analysis
data=ROImatrix';
[num_record,data_len]=size(data);
amplitudes = cell(num_record,1);

for j=1:num_record
    spike_data = data(j,:);
    sigma = std(data(j,:)); 
    threshold = mean(data(j,:)) + sigma_threshold * sigma;
    cd('../');
    spike_index = spike_detection(spike_data, threshold);
    cd(savename);
    amplitude = data(j, spike_index);
    timestamp{j} = sparse(j,spike_index,amplitude);
    amplitudes{j} = amplitude;
end

figure;
hold on;
for i = 1:num_record
    if length(timestamp{i} ~= 0)
        plot((num_record - i) + timestamp{i}(i,:), '-', 'linewidth', 1,'Color',colors(i,:));
    end
end
set(gca, 'YTick', [], 'ycolor', 'k');
ylabel('All spikes','Fontname','Arial','FontSize', 10);
xlabel('Frame','Fontname','Arial','FontSize', 10);
hold off;
print('-dtiff',strcat(savename, '_allspikes.png'));

mid=round(data_len/2);
for m = 1:num_record
    min1 = min(data(m,1:100));
    min2 = min(data(m,(mid-50):(mid+50)));
    min3 = min(data(m,(data_len-100):(data_len)));
    mins = [min1 min2 min3];
    F0(m,1) = mean(mins);
    amplitudes_dff{m} = (amplitudes{m} - F0(m,1))./F0(m,1);
    deltaFtoF0(m,1) = mean(amplitudes_dff{m});  
end
ave_dff = mean(deltaFtoF0);
std_dff = std(deltaFtoF0);

for k = 1:num_record
    [~, num_spikes{k}] = size(amplitudes{k});
    neuron_freq(k,1) = num_spikes{k} / (filetime/60);
end
ave_freq = mean(neuron_freq);
std_freq = std(neuron_freq);

figure;
reflclr=[0.8,0.8,0.8];
refclr='k';
F0plot=F0';
for i=1:num_record
    dFF0(:,i) = (ROImatrix(:,i)-F0plot(:,i))./F0plot(:,i);
end
b=1:data_len;
t=b*timeinterval;
plot(t,dFF0,'Color',reflclr,'LineWidth',0.5,'LineStyle','-');hold on;
avg=mean(dFF0');
plot(t,avg,'-o','Color',refclr,'LineWidth',2,'MarkerSize',2); 
ylim([-0.5 1.5]);
xlim([0 filetime]);
ylabel('Normalized GCaMP Intensity', 'Fontname','Arial', 'FontSize', 10);
xlabel('Time (s)', 'Fontname','Arial','FontSize', 10);
print('-dtiff',strcat(savename, '_normalizedplot.png'));

header = cell(1,3);
for i=1:num_record
    header{1,i}=['Neuron ' num2str(i)];
end

writecell([header;num2cell(ROImatrix)],strcat(savename, '_data.xlsx'),'Sheet','GCaMP Intensity');
writecell([header;num2cell(dFF0)],strcat(savename, '_data.xlsx'),'Sheet','Normalized Intensity');

cd('../');

spectralanalysis = {'Frequency', 'Period', 'Power', 'Percentage'};

Summary = cell(2+num_record,4);
Summary{1,2}='Frequency (spikes/min)';
Summary{1,3}='Average amplitude';
Summary{1,4}='Std dev amplitude';
Summary{2,1}='Synchronous';
Summary{2,2}=Sync_rate;
Summary{2,3}=mean_Sync_amplitude_dFF;
Summary{2,4}=std_Sync_dff;
spikes = cell2mat(num_spikes);
for i=1:num_record
    [peaks, ~] = fft_neuron(ROImatrix(:,i),filename,timeinterval);
    spectralanalysis = vertcat(spectralanalysis,num2cell(peaks));
    
    Summary{2+i,1}=['Neuron ' num2str(i)];
    Summary{2+i,2}=double(spikes(i))/(filetime/60);
    amps = cell2mat(amplitudes_dff(i));
    Summary{2+i,3}=mean(amps);
    Summary{2+i,4}=std(amps);
end

cd(savename);
writecell(spectralanalysis,strcat(savename, '_data.xlsx'),'Sheet','Spectral Analysis');
writecell(Summary,strcat(savename, '_data.xlsx'),'Sheet','Statistics Summary');
cd('../');
end

function [spikes_index] = spike_detection(spike_data,threshold)
spike_index = find((spike_data)>threshold);
spike_index_array = zeros(1,length(spike_index));
i=1;

while ( i < length(spike_index))
    n=1;
    temp_var = 0;
    while((spike_index(i) == spike_index(i+n)-n) && i+n<length(spike_index))
        n=n+1;
        if(temp_var<spike_data(spike_index(i+n)))
            temp_var = spike_data(spike_index(i+n));    
            spike_index_array(1,i) = spike_index(i+n);
        end 
    end 
    i=i+n;
end

spikes_index = spike_index_array(spike_index_array~=0);
end 

