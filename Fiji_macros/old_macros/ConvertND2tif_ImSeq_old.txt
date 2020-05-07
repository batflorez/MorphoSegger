// This macro runs Bio-formats to convert ND2 files to tif and perform operations 
// to prepare data for SuperSegger and Morphometrics analysis
// Author: Andres Florez     April 27, 2020

inDir = getDirectory("Choose input Directory");
outDir = getDirectory("Choose output Directory");
setBatchMode(true);
//print("-- Done --");
filenames = getFileList(inDir);
pattern = ".*"; // for selecting all the files in the folder
timePrefix="_t";
extension=".nd2"
// different examples of regex-based selection 
//pattern = "01-03.*Pos [3-9].*";
//pattern = ".*Pos [7-9].*";
//pattern = "01-02.*";

//Add additional channels if necessary, always test on macro recorder

count = 0;
for (i = 0; i < filenames.length; i++) {
	currFile = inDir+filenames[i];
	if(endsWith(currFile, extension) && matches(filenames[i], pattern)) { // process ND2 matching regex, change here for other file types
		//open(currFile);
		count++;
        run("Bio-Formats Importer", "open=currFile autoscale color_mode=Default rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
        fluor = getImageID();
        selectImage(fluor);
        fluorName= getTitle();
        titleFluor = replace(fluorName, extension, "xy");
        rename(titleFluor);
        run("Gaussian Blur...", "sigma=1 stack");
        run("Image Sequence... ", "format=TIFF name=["+titleFluor+timePrefix+"] save=outDir");
        //saveAs("Tiff", outDir+titleFluor);
        close();
        phase = getImageID();
        selectImage(phase);
        phaseName= getTitle();
        titlePhase = replace(phaseName, extension, "xy");
        rename(titlePhase);
        run("Image Sequence... ", "format=TIFF name=["+titlePhase+timePrefix+"] save=outDir");
        close();
        
        //title = replace(filenames[i], ".nd2", "");
        //rename(title);
        //selectWindow(filenames[i]+" - C=0");
        //run("Image Sequence... ", "format=TIFF name=glucose01 use save=outDir");
        //close();
        //selectWindow(title+" - C=1");
        //run("Gaussian Blur...", "sigma=1 stack");
        //saveAs("Tiff", outDir+title);
		//close(); // close C0
	}
}
print("Number of files processed: "+count);

