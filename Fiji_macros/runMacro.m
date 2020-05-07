% prepareFiles: This function runs MIJ to open an instance of Fiji without
% graphical interface and runs a macro without arguments.
% Andres Florez - 04/29/2020
% andrewflorez@gmail.com
% Harvard University

function runMacro(macroDir,args)

% args refers mostly to file paths to prevent the macro asking for a
% directory.

Miji(false); %open MIJ withouth Fiji's graphical interface
IJ=ij.IJ(); %Create an instance of Fiji's methods

IJ.runMacroFile(java.lang.String(macroDir),java.lang.String(args));   %run macro

MIJ.exit() %closes MIJ

end