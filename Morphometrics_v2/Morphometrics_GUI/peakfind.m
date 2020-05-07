function varargout =peakfind(y_data,varargin)
%
%PEAKFIND general 1D peak finding algorithm
%Tristan Ursell, 2011.
%
%peakfind(y_data)
%peakfind(x_data,y_data)
%peakfind(x_data,y_data,upsam)
%peakfind(x_data,y_data,upsam,gsize,gstd)
%peakfind(x_data,y_data,upsam,htcut,'cuttype')
%peakfind(x_data,y_data,upsam,gsize,gstd,htcut,'cuttype')
%
%[peakspos]=peakfind()
%[xout,yout,peakspos]=peakfind()
%
%This function finds peaks without taking first or second derivatives,
%rather it uses local slope features in a given data set.  The function has
%four basic modes.
%
%   Mode 1: peakfind(x_data,y_data) simply finds all peaks in the data
%   given by 'xdata' and 'ydata'.
%
%   Mode 2: peakfind(x_data,y_data,upsam) finds peaks after up-sampling
%   the data by the integer factor 'upsam' -- this allows for higher
%   resolution peak finding.  The interpolation uses a cubic spline that
%   does not introduce fictitious peaks.
%
%   Mode 3: peakfind(x_data,y_data,upsam,gsize,gstd) up-samples and then
%   convolves the data with a Gaussian point spread vector of length
%   gsize (>=3) and standard deviation gstd (>0).  The convolution option 
%   only works with upsam>=2.
%
%   Mode 4: peakfind(x_data,y_data,upsam,htcut,'cuttype') up-samples the
%   data (upsam=1 analyzes the data unmodified).  The string 'cuttype' can
%   either be 'abs' (absolute) or 'rel' (relative), which specifies a peak 
%   height cutoff which is either:
%
%       'abs' - finds peaks that lie an absolute amount 'htcut' above the
%       lowest value in the output y_data set.
%
%            for (htcut > 0) peaks are found if
%
%            peakheights - min(yout) > htcut
%       
%       'rel' - finds peaks that are an amount 'htcut' above the lowest
%       value in the data set, relative to the full change in y-input
%       values.
%
%           for (0 < htcut < 1) peaks are found if
%
%           (peakheights-min(yout))/(max(yout)-min(yout)) > htcut
%
%Upsampling and convolution allows one to find significant peaks in noisy
%data with sub-pixel resolution.  The algorithm also finds peaks in data 
%where the peak is surrounded by zero first derivatives, i.e. the peak is 
%actually a large plateau.  
%
%The function outputs the x-position of the peaks in 'xpeaks' or the
%processed input data in 'xout' and 'yout' with 'peakspos' as the indices
%of the peaks, i.e. xpeaks = xout(peakspos).  
%
%If you want the algorithm to find the position of minima, simply input 
%'-y_data'.  Peaks within half the convolution box size of the boundary 
%will be ignored (to avoid this, pad the data before processing).
%
% Example:
%
% x_data = -50:50;
% y_data =(sin(x_data)+0.000001)./(x_data+0.000001)+1+0.025*(2*rand(1,length(x_data))-1);
%
% [xout,yout,peakspos]=peakfind(x_data,y_data,4,6,2,0.2,'rel');
%
% plot(x_data,y_data,'r','linewidth',2)
% hold on
% plot(xout,yout,'b','linewidth',2)
% plot(xout(peakspos),yout(peakspos),'g.','Markersize',30)
% xlabel('x')
% ylabel('y')
% title(['Found ' num2str(length(peakspos)) ' peaks.'])
% box on

if and(nargout~=1,nargout~=3)
    error('Incorrect number of output arguments.')
end

%*****************************
%*** PARSE VARARGIN **********
%*****************************

%parse varargin input
if any([nargin==4,nargin==6,nargin>7])
    error('Incorrect number of input arguments.  See "help peakfind".')
end

if nargin==1
    x_data=1:numel(y_data);
else
    if nargin>1
        %peakfind(x_data,y_data)
        x_data=y_data;
        y_data=varargin{1};
    end
        
    %check if input data are same length
    if length(x_data)~=length(y_data)
        error('Input data vectors must be the same length.')
    end
    
    %check to make sure data is all real
    if any([~isreal(x_data),~isreal(y_data)])
        error('All input data must be real valued.')
    end
    
    %set the second varargin to 'upsam'
    if nargin>2
        upsam=round(varargin{2});
        
        if upsam<1
            error('Upsampling factor must be an integer greater than or equal to 1.')
        end
    else
        upsam=0;
    end
    
    if nargin==5
        if ischar(varargin{4})
            %peakfind(x_data,y_data,upsam,htcut,'cuttype')
            htcut=varargin{3};
            cuttype=varargin{4};
        else
            %peakfind(x_data,y_data,upsam,gsize,gstd)
            gsize=varargin{3};
            gstd=varargin{4};
        end
    end
    
    %peakfind(x_data,y_data,upsam,gsize,gstd,htcut,'cuttype')
    if nargin==7
        gsize=varargin{3};
        gstd=varargin{4};
        htcut=varargin{5};
        cuttype=varargin{6};
    end
end


%**************************
%*** CHECK INPUT **********
%**************************

%check if input data is monotonic or flat (no peaks)
if or(all(diff(y_data)<=0),all(diff(y_data)>=0))
    %disp('Warning: Input data is monotonic.')
    if nargout==1
        %[peakspos]=peakfind()
        varargout(1)={[]};
    elseif nargout==3
        %[xout,yout,peakspos]=peakfind()
        varargout(3)={[]};
        if nargin==1
            varargout(2)={y_data};
            varargout(1)={1:numel(y_data)};
        else
            varargout(2)={varargin{1}};
            varargout(1)={y_data};
        end
    else
        error('Incorrect number of output arguments.')
    end   
    return
end

%check input values
if exist('gsize','var')
    if gsize<3
        error('Gaussian window must have at least 3 points.')
    end
    
    if mod(gsize,2)==0
        gsize=gsize+1;
    end
end

if exist('gstd','var')
    if gstd<=0
        error('Gaussian STD must be greater than zero.')
    end
end

if exist('htcut','var')
    if htcut<=0
        error('The variable ''htcut'' must be greater than zero.')
    end
end

if exist('cuttype','var')
    if and(~strcmp('rel',cuttype),~strcmp('abs',cuttype))
        error('Invalid string entry in function call for variable ''cuttype''.')
    end
end

%upsample data
if upsam>1
    xup_data=linspace(x_data(1),x_data(end),upsam*length(x_data));
    yup_data=interp1(x_data,y_data,xup_data,'pchip');
else
    xup_data=x_data;
    yup_data=y_data;
end

%create convolution filter
if and(exist('gsize','var'),exist('gstd','var'))
    %define convolution filter
    filt1 = exp(-linspace(-gsize/2,gsize/2,gsize).^2/(2*gstd^2));
    conv1 = filt1/sum(filt1);
    norm1 = conv(ones(1,numel(yup_data)),conv1,'same');
    %convolve peak data with small STD gaussian
    ycv_data=conv(yup_data,conv1,'same')./norm1;
else
    ycv_data=yup_data;
end

xcv_data=xup_data;

%construct peak finding vectors
lpts=ycv_data(1:end-2);
mpts=ycv_data(2:end-1);
rpts=ycv_data(3:end);

%number of upsampled data points
tipn=length(ycv_data);

%find peaks
Rpeaks=zeros(1,tipn);
Lpeaks=zeros(1,tipn);
minpeaks=zeros(1,tipn);
fpeaks=zeros(1,tipn);

Rpeaks(2:end-1)=(mpts>=lpts).*(rpts<mpts); %triggers just right of peak
Lpeaks(2:end-1)=(mpts>lpts).*(rpts<=mpts); %triggers just left of peak
minpeaks(2:end-1)=(mpts<=lpts).*(rpts>mpts);

%use minimums and multiple maximums to find best guess of peak
minlist=find(minpeaks==1);

%triggers when there is a flat landscape with a single trough
if and(all(~Rpeaks),all(~Lpeaks))
    warning('Found only a trough.')
    npeaks=0;
    peakspos=[];
    xcv_data=x_data;
    ycv_data=y_data;
    return
end

%find and average doubly-occupied peaks
for j=1:length(minlist)+1
    %find bounds between which to look for two pseudo-peaks
    if j==1
        lbound=1;
        if isempty(minlist)
            rbound=tipn;
        else
            rbound=minlist(j);
        end
    elseif j==length(minlist)+1
        lbound=minlist(j-1);
        rbound=length(ycv_data);
    else
        lbound=minlist(j-1);
        rbound=minlist(j);
    end
    
    %look for peaks
    peakmask=zeros(1,tipn);
    peakmask(lbound:rbound)=1;
    
    Lpeak=find(Lpeaks.*peakmask==1,1);
    Rpeak=find(Rpeaks.*peakmask==1,1);
    fpeaks(round((Lpeak+Rpeak)/2))=1;
end

%do a final check against the peaks in the original up-sampled data
peakspos=find(fpeaks==1);
npeaks=sum(fpeaks);

%remove peaks if they are not high enough
if and(npeaks>1,exist('cuttype','var'))
    ymin=min(ycv_data);
    ymax=max(ycv_data);
    
    if strcmp('rel',cuttype)
        cutpeaks=((ycv_data(peakspos)-ymin)/(ymax-ymin))<htcut;
    else
        cutpeaks=ycv_data(peakspos)<(htcut+ymin);
    end
    npeaks=npeaks-sum(cutpeaks);
    peakspos=peakspos(cutpeaks==0);
end

%format output arguments
if nargout==1
    varargout(1)={xcv_data(peakspos)};
else
    varargout(3)={peakspos};
    varargout(2)={ycv_data};
    varargout(1)={xcv_data};
end

