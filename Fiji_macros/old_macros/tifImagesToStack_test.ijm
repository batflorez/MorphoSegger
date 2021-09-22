
// This macro converts  all tif image sequences in a folder to stack.
// The macro searches for xy folders and process fluor1, phase and seg subfolders.
// Author: Andres Florez     April 29, 2020. Harvard University


//Select input
//inDir = getArgument(); // get path of Analysis folder from Matlab as an argument

inDir = '/Users/aflorez/Dropbox/Projects_postdoc17/Manuscript_1/Data/Calarco/bAF217_S750gluc_3min_061421_doubleCheckThisExperiment/Analysis_c/'
dirnames = getFileList(inDir); 

//Batch Mode
setBatchMode(true);

count1 = 0;
count2 = 0;

for (n = 0; n < dirnames.length; n++) {    // lists all the directories under Analysis folder
        currDir1 = inDir+dirnames[n];
        list= getFileList(currDir1);
        for (k=0; k<list.length; k++) {
            currDir2 = inDir+dirnames[n]+list[k];  // lists the directories under xy folders
            count1++;
            print(currDir2);
            if (endsWith(list[k],"fluor1/")){
                run("Image Sequence...", "dir=currDir2 sort");
                //name1=getTitle();
                //saveAs("Tiff", currDir2+name1+"_xy"+(n)+".tif"); // save stack as .tif
                saveAs("Tiff", currDir2+"_xy"+(n)+".tif"); // save stack as .tif
                close();
            }

            //if (endsWith(list[k],"phase/")){
            //    run("Image Sequence...", "dir=currDir2 sort");  // save stack as .tif
                //name2=getTitle();
                //saveAs("Tiff", currDir2+name2+"_xy"+(n)+".tif");
            //    saveAs("Tiff", currDir2+"_xy"+(n)+".tif");
            //    close();
            //}
            //
            //list2= getFileList(currDir2);  // delete all the single image files (.tif)
            //Array.print(list2);
            //for (c=0; c<list2.length; c++) {
            //    count2++;
            //    if (startsWith(list2[c],"t") && endsWith(list2[c], ".tif")){  //this only works if all the images always start with time variable (t) and ends with .tif
            //        ok=File.delete(currDir2+list2[c]);
            //    }
            } 
         }   
//}
count3 = count2-count1;
print("Number of folders processed "+count1);
print("Number of image files deleted: "+count3);


