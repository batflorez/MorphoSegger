function params = loadParams( paramsName )
% loadConstants : loads the parameters for the superSegger/trackOpti package.
%This function was adapted from SuperSegger - Andres Florez 04/26/20

% gets the list of all possible constants in the Parameter_files folder
[possibleConstants,list, filepath] = getParamsList();



indexConst = find(strcmpi({possibleConstants.name},[paramsName ,'.mat']));
if ~isempty(indexConst)
     constFilename = possibleConstants(indexConst).name;
     params = load([filepath,constFilename],'f_*','v_*');
else
    disp('loadParams: Parameters not found. Aborting. ');
    disp(['Possible constants']);
    disp(list');
    params = [];
    return;
end

