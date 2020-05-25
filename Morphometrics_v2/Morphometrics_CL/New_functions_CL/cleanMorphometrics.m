function cleanMorphometrics(dirname)

% This function deletes Gparent and Gsegt (temporary) files when
% Morphometrics stops due to errors.

% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University


dirname =fixDir(dirname);           % The analysis folder
contents = dir([dirname,'xy*']);    % List all xy folders
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

for k = 1:num_xy
    delete([dirnamelist{k},filesep,'*_Gsegt.tif']);
    delete([dirnamelist{k},filesep,'*_Gparent.tif']);    
end
clearvars

