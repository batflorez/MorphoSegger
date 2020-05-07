% runMorphometrics                              %
% Input: Analysis folder, parameter file 
%
% This function runs morphometrics_mask_cl on multiple xy points in parallel
% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University

function run_parallel(dirname, params)

%dirname=pwd;
%dirname=fixDir(dirname);
contents = dir([dirname,'xy*']); %List al xy folders
num_dir_tmp = numel(contents);
nxy = [];
num_xy = 0;

%Creates a struct for all seg directories 
dirnamelist=cell(1,num_dir_tmp);
for i = 1:num_dir_tmp
    if (contents(i).isdir) && (numel(contents(i).name) > 2)
    num_xy = num_xy+1;
    nxy = [nxy, str2double(contents(i).name(3:end))];
    dirnamelist{i} = [dirname,contents(i).name,filesep,'seg'];
    end
end

%List seg mask files and runs morphometrics in parallel
parfor (k = 1:num_xy)
    mask = dir([dirnamelist{k},filesep,'*.tif']);
    maskFile_tmp= [mask.folder,filesep, mask.name];
    morphometrics_mask_cl(maskFile_tmp,params,[],1); 
end

% move morphometrics files to independent folder
for q = 1:num_xy
    movefile([dirnamelist{q},filesep,'*.mat'], [dirname,contents(q).name,filesep,'morphometrics']);
end

disp('Parallel Morphometrics done.')
% %workers=6;
% if  workers % shutting down parallel pool     
%     delete(poolobj);
% end

end