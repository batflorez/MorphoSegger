// This macro converts  all tif images in a folder to stack
// Author: Andres Florez     April 29, 2020

//run("Image Sequence...", "sort");

// Make a loop through the fluor and phase directories
// run the image Sequence thingy

//Select input, and create output directory
inDir = getDirectory("Choose input Directory");
dirnames = getFileList(inDir); //list of files

//Batch Mode
setBatchMode(true);

pattern = "xy[0-9]"; // selection xy folders
// different examples of regex-based selection 
//pattern = "01-03.*Pos [3-9].*";

count = 0;
for (i = 0; i < dirnames.length; i++) {
	//if(matches(dirnames[i], pattern)) { // process matching regex, 
        currDir = inDir+dirnames[i];
        count++;
        print(currDir);
        list= getFileList(currDir);
        for (k=0; k<list.length; k++) {
            currDir2 = inDir+dirnames[i]+list[k];
            print(currDir2);
            if (endsWith(list[k],"fluor1/"))
            print(currDir2);
            
         }
        //print(tifDir);
        //processFile(currDir);
        //Add additional channels if necessary, always test on macro recorder
        //saveAs("Tiff", outDir+title);
		//close(); 
	//}
}


// I need another for to go through xy1
// then do phase and fluor1
print("Number of folders processed: "+count);

function processFile(path) {
    if (endsWith(path, ".tif")) {
        open(path);
        run("Image Sequence...", "sort");
        save(path);
        close();
   }
}



