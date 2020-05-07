%Tristan Ursell
%April 2012
%Function to remove false positive regions from segmented bacterial images
%
%Lin image must be bwlabel or segmented BW image.
%Iin image is the grayscale from which Lin was generated.
%
%typeq = 'fluor_i','flour_p',or 'phase'.
%

function Lout=falsepos(Lin,Iin,int_rej)

global f_hmin

%check intensities of regions between border and interior
L1=bwlabel(Lin);
L1_perim=bwperim(L1);

%find the list of regions that appear in both interior and perimeter forms
%this is to handle the very rare case that a region only has perimeter
%pixels
L1_list=L1-L1_perim.*L1;
v1=unique(L1(L1>0));
v2=unique(L1_list(L1_list>0));
setdiff1=setdiff(v1,v2);

for i=1:length(setdiff1)
    L1(L1==setdiff1(i))=0;
end

%reform inputs
L1=bwlabel(L1);
L1_perim=bwperim(L1);
Lout=L1;

%find statistics for segmented regions
temp1=regionprops(~L1_perim.*L1,Iin,'MeanIntensity','PixelValues');
temp2=regionprops(L1_perim.*L1,Iin,'MeanIntensity');

interiormean=[temp1.MeanIntensity]; %segmented interior pixel means
perimmean=[temp2.MeanIntensity]; %segemented border pixel means

interiorstd=zeros(1,length(temp1)); %segmented interior pixel STDs
for j=1:max(L1(:))
    interiorstd(j)=std(temp1(j).PixelValues);
end

%calculate border / intertior ratios (should be greater than 1)
brd_ratio=perimmean./interiormean;
%{
if strcmp(typeq,'fluor_i')
    brd_ratio=interiormean./perimmean;
else
    brd_ratio=perimmean./interiormean;
end
%}

%nois / mean ratio
noise_ratio=interiorstd./interiormean;

%remove false positive regions based on border / interior ratios
cond1=brd_ratio<=(1+int_rej*f_hmin); %if the mean of the region perimeter pixels is inappropriate compared to the mean of the interior pixels
cond2=noise_ratio<f_hmin/2; %if the STD of the interior pixel regions is too small
false_pos=find(or(cond1,cond2));
for j=1:length(false_pos)
    Lout(L1==false_pos(j))=0;
end

%turn to logical, and remove spurs
Lout=bwmorph(Lout>0,'spur',Inf);
















