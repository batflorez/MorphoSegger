%Image Stack Saver
%Tristan Ursell, February 2012
%
% image_save(Im1,basename)
% image_save(Im1,basename,fmax)
%
%For use with appended image stacks, e.g. most commonly stacked TIF files.
%Some operating systems throw an error indicating that the file is
%inaccessible at the time of write, which leads to a code fault and can be
%very frustrating.  This simple script fixes that issues, mainly on the
%Windows operating system.
%
% Im1 = the current image (matrix) you want to add to the stack
% 
% basename = is a string that specifices the file name, and potentially the
% path of the stack to save to.  Best practice is the put '.tif' at the end
% of the file name.  
%
% fmax = maximum number of allowed failures, i.e. if the script attempts to
% write to the file 'fmax' times and fails, it will give up.  This prevents
% entering an infinite loop.  The default value is fmax = 10.
%
% Any of the 'imwrite' parameters can be modified below on line 38.
%

function image_save(Im1,basename,varargin)

if ~isempty(varargin)
    fmax=varargin{1};
else
    fmax=10;
end

er1=0;
f=0;
while and(er1==0,f<fmax)
    try
        imwrite(Im1,basename,'Compression','none','WriteMode','append');
        er1=1;
    catch
        f=f+1;
        pause(0.1*rand)
    end
end

if f==fmax
    error(['File: ' basename ', remained inaccessible after ' num2str(fmax) ' attempts.'])
end


