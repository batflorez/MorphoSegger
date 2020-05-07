This is the Matlab (c) script version of Morphometrics.


************************ NEW VERSION 2.0 Apr 2020 **************************************

This revised version of Morphometrics (from Tristan Ursell) works with Matlab 2019b 

It works now with Mac and Linux (all bugs have been fixed). To make it work for each distribution, unzip the MM_mexfiles for the appropriate OS and add them to the path. (Delete the other OS mexfiles from the path to prevent problems). In the MAC version, make sure to allow the mex files be opened on privacy settings.

simply_segment_cl.m had several bug fixes. The false pos function used for the command line version was the the falsepos.m file instead of the one written on simply_segment_cl.m

It saves automatically the cells_table_ID in the command line version, previously only in the GUI.

morphometrics_mask_cl script for binary masks. It can run in parallel (multiple cores) for different XY positions when using the run_parallel.m script. It requires parallel computing toolbox.

Several bugs have been fixed in view_contours function and scripts for command line processing (cl scripts).

v2struct function was added to easily access and change parameter values as a structure when running a pipeline.

This new version can run in a Pipeline with SuperSegger Segmentation.

To start, type morphometrics to initialize the GUI (with all files on path) 


Andres Florez
Harvard University
andrewflorez@gmail.com

***************************PREVIOUS VERSION ****************************************

You need Matlab 2009 or higher installed, including the Mapping and Image Processing Toolboxes.

Generally, this program works better on PCs than Macs.

Keep all these files together.

To launch the program, simply navigate Matlab to the directory containing these files, and then type 'morphometrics' in the command window.

Alternatively, add the folder path by going to "File" -->  "Set Path..."  --> "Add with subfolders..."  -->  choose the "Morphometrics" directory --> "Save".

**Please read the Help file (by pressing 'Help' in the program).**

Thanks.
Tristan
tsu@uoregon.edu

