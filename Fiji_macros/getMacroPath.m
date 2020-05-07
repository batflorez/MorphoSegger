
function filepath = getMacroPath()

FulllocationOfFile = mfilename('fullpath');
fileSepPosition = find(FulllocationOfFile==filesep,1,'last');
filepath = FulllocationOfFile ( 1 :fileSepPosition-1);
filepath = [filepath,filesep];

end

