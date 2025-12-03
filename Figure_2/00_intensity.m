function [] = Intensity(filename,filetime,diameter,sigma_threshold)
stack_info = imfinfo(filename);

[~,savename,~]=fileparts(filename);
mkdir(savename);

image_width = stack_info(1).Width;
image_height = stack_info(1).Height;
image_count = length(stack_info);
image_stack = zeros(image_height,image_width,image_count,'uint16');
timeinterval = filetime/image_count;

[t]=(0:image_count);
v=zeros(image_count,1);
for i=1:image_count
   v(i)=t(i)*timeinterval;
   frame = imread(filename,'Index',i);
   image_stack(:,:,i)=frame(:,:,2);
end

cd(savename);

max_image=uint16(max(double(image_stack),[],3));
screen_size = get(0,'ScreenSize');
figure('Position',[screen_size(3)*0.2 0 image_width image_height]);
set(gca,'Position',[0 0 1 1]);
imshow(imadjust(max_image,stretchlim(max_image)));
[col,row,~]=impixel;

colors=parula(length(col));

ROI_intensity=zeros(image_count,length(row));
for k=1:length(row)
   ROI_intensity(:,k)=mean(mean(image_stack(max(1,row(k)-diameter):min(row(k)+diameter,image_height),max(1,col(k)-diameter):min(col(k)+diameter,image_width),:)));
end

screen_res=get(0,'ScreenPixelsPerInch');
figure('Position',[screen_size(3)*0.2 0 image_width image_height],'PaperPosition',[0.25 2.5 image_width/screen_res image_height/screen_res],'PaperUnits','inches');
set(gca,'Position',[0 0 1 1]);
imshow(imadjust(max_image,stretchlim(max_image)));
axis off;
set(gca,'XTick',nan,'YTick',nan);
hold on;
for k=1:length(row)
   rectangle('Position',[col(k)-diameter,row(k)-diameter,2*diameter+1,2*diameter+1],'Curvature',[1 1],'EdgeColor',colors(k,:),'LineWidth',2);
   hold on;
end
hold off;
print('-dtiff',strcat(savename, '_selections.png'));

save(strcat(savename, '_selections.mat'),'ROI_intensity','col','row','diameter','filename','timeinterval','colors');

b=1:length(ROI_intensity(:,1));
t=b*timeinterval;
figure;
for k=1:length(row)
    plot(v,ROI_intensity(:,k),'LineWidth',0.5,'Color',colors(k,:));
    hold on;
end
hold on;
avg=mean(ROI_intensity');
plot(t,avg,'Color','k','LineWidth',1.5);
xlim([0 timeinterval*image_count])
ylabel('GCaMP Intensity','Fontname','Arial','FontSize', 10);
xlabel('Time (s)', 'Fontname','Arial','FontSize', 10);
hold off;
print('-dtiff',strcat(savename, '_plot.png'));

cd('../');
IntensityAnalysis(filename,ROI_intensity,filetime,sigma_threshold,colors);
close all;
end
