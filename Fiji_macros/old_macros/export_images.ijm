/*
 * Complex Format EXPORT MACRO
 * By Olivier Burri @ EPFL - SV - PTECH - BIOP
 * Given a folder, extracts all series inside all multi-file files with given extension in new folders 
 * Last edit: 13.02.2017
 */


 ////////////////////// SET PARAMETERS //////////////////////
 ////////////////////////////////////////////////////////////
 
 
// Set the extension you would like this macro to work with.
// Do not add a . at the beggining
extension = "nd2";  //eg "lif", "vsi", etc...


// Set to true if you want all planes of the image to be saved individually
// if set to false, it will save each series as a stack. 
is_save_individual_planes = true; // set to either true or false

// Padding for the naming of the series if you want to save all 
// Images individually
pad = 3; // 0 means no padding. 2 means '01', '02' etc...


 //////////////////// END SET PARAMETERS ////////////////////
 ////////////////////////////////////////////////////////////
 




// Beggining of macro. You should now have anything to edit after this line. 

dir = getDirectory("Select a directory containing one or several ."+extension+" files.");

files = getFileList(dir);


setBatchMode(true);
k=0;
n=0;

run("Bio-Formats Macro Extensions");
for(f=0; f<files.length; f++) {
	if(endsWith(files[f], "."+extension)) {
		k++;
		id = dir+files[f];
		Ext.setId(id);
		Ext.getSeriesCount(seriesCount);
		print(seriesCount+" series in "+id);
		n+=seriesCount;
		for (i=0; i<seriesCount; i++) {
			run("Bio-Formats Importer", "open=["+id+"] color_mode=Default split_channels view=Hyperstack stack_order=XYCZT series_"+(i+1));
			fullName	= getTitle();
			dirName 	= substring(fullName, 0,lastIndexOf(fullName, "."+extension));
			fileName 	= substring(fullName, lastIndexOf(fullName, " - ")+3, lengthOf(fullName));
			File.makeDirectory(dir+File.separator+dirName+File.separator);

			print("Saving "+fileName+" under "+dir+File.separator+dirName);
			
			getDimensions(x,y,c,z,t);
			
			
			if(is_save_individual_planes) {
				save_string = getSaveString(pad);
				print(dir+File.separator+dirName+File.separator+fileName+"_"+save_string+".tif");
                //run("Image Sequence... ", "format=TIFF name=["+fileName+"] digits="+pad+" save=["+dir+File.separator+dirName+File.separator+"]");
                run("Image Sequence... ", "format=TIFF name=["+fileName+"] digits="+pad+" save=["+dir+dirName+File.separator+"]");
				
			} else {
                //saveAs("tiff", dir+File.separator+dirName+File.separator+fileName+"_"+(i+1)+".tif");
                saveAs("tiff", dir+dirName+File.separator+fileName+"_"+(i+1)+".tif");
			}
			run("Close All");
		}
	}
}
Ext.close();
setBatchMode(false);
showMessage("Done with "+k+" files and "+n+" series!");


function getSaveString(pad) {
	str ="";
	getDimensions(x,y,c,z,t);
	if(t > 1)  str+="t_"+IJ.pad(1,pad);
	if(z > 1)  str+="z_"+IJ.pad(1,pad);
	if(c > 1)  str+="c_"+IJ.pad(1,pad);
	
	return str;
}