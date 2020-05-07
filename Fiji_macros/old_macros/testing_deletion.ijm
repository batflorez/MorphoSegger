// This macro converts  all tif images in a folder to stack
// The macro enters xy folders and searches for fluor1 and phase subfolders
// Author: Andres Florez     April 29, 2020


//Select input
inDir = getDirectory("Choose input Directory");
dirnames = getFileList(inDir); //list of files

//Batch Mode
setBatchMode(true);

count = 0;
for (i = 0; i < dirnames.length; i++) {    // lists all the directories under Analysis folder
        currDir1 = inDir+dirnames[i];
        //print(currDir);
        list= getFileList(currDir1);
        for (k=0; k<list.length; k++) {
            currDir2 = inDir+dirnames[i]+list[k];  // lists the directories under xy folders
            //print(currDir2);
            //save fluorescence images as stack
            list2= getFileList(currDir2);
            for (c=0; c<list2.length; c++) {
                if (endsWith(list2[c],".tif")){
                    ok=File.delete(currDir2+list2[c]);
                }
            }
        }
    
}

print("Number of folders processed: "+count);
//Array.print(list2); //
//one trick i can use is i can save the stack as .tiff and delete the other files.
//I need to test deleting files so i can start wi