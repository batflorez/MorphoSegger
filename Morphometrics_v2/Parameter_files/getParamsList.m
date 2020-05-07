function [possibleConstants,  list, filepath] = getParamsList()
% getConstantsList : gets all constants files from the settings directory.
%This function was adapted from SuperSegger - Andres Florez 04/26/20

FulllocationOfFile = mfilename('fullpath');
fileSepPosition = find(FulllocationOfFile==filesep,1,'last');
filepath = FulllocationOfFile ( 1 :fileSepPosition-1);
filepath = [filepath,filesep];
possibleConstants = dir([filepath,'*.mat']);

list = {};

for i = 1 : numel (possibleConstants)
    cName = possibleConstants (i).name;
    possibleConstants(i).resFlag =cName (1:end-4);
    list{i} = possibleConstants(i).resFlag;
end

end