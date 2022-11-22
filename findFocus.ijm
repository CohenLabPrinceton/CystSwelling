
directory = "//latrobe/CohenLabArchive/Gawoon & Isaac/cyst images/060122_ctrl/CFTR/"

//get information from the first .tif file about the dimensions
open(toString(directory+"/cropped/cyst1.tif")); 
Stack.getDimensions(width, height, channels, slices, frames);

firstTP=1; //first t point
lastTP=frames; //last t point
stackSize = slices; //how many timepoints are there in the stack?

close();

numCysts= 72; // number of cysts in this experiment folder

// go through each z-stack file and find the focused slide index using "Find Focused Slices"
// for each cyst, the index of focused image at each timepoint is saved in a .csv file 
for (nc=1; nc<=numCysts; nc++) { 

	
	cystnum = toString(nc);
	open(toString(directory+"cropped/cyst"+toString(nc)+".tif")); 

	cystname = "cyst"+cystnum+".tif";
	slices = newArray(lastTP);
	selectWindow(cystname);
	
	
	for (i=1; i<=lastTP; i++) {
	
		newFrame = toString(i);
		selectWindow(cystname);
		run("Duplicate...", "duplicate frames="+newFrame);
		run("Find focused slices", "select=100 variance=0.8 edge log");
		
		slice = getResult('Slice',0); // Find Focused Slides finds the slice in focus
		slices[i-1] = slice + 1; 
	
		
		selectWindow("cyst"+cystnum+"-1.tif");
		close();
		selectWindow("Focused slices of cyst"+cystnum+"-1.tif_100.0%"); //
		close();
	
	}

	
	for (i=1; i<=lastTP; i++) {
		
		setResult("Focus", i-1, slices[i-1]); // Save the index of the focused slide
	}
	
	updateResults();
	saveAs("Results", directory+"cyst"+cystnum+".csv"); 
	
	// close all files
	if (isOpen("Results")) {
         selectWindow("Results"); 
         run("Close" );
	}
    if (isOpen("Log")) {
         selectWindow("Log");
         run("Close" );
    }
     while (nImages()>0) {
          selectImage(nImages());  
          run("Close");
    }
}

// use the .csv file from the previous step to find, extract, and concatenate the in-focus slides for

File.makeDirectory(directory+File.separator+"focused_cysts");

for (nc=1; nc<=numCysts; nc++) {

	cystnum = toString(nc);
	
	open(toString(directory+"cropped/cyst"+toString(nc)+".tif")); 
	sliceInfo = "cyst"+cystnum+".csv";
	open(directory+sliceInfo);
	
	for (i=1; i<=lastTP; i++) {
		
		stack = getResult('Focus',i-1); // Read the index of the focused slide
		
		if (stack == stackSize + 1){
			stack = stackSize;
		}

		
		zframe = toString(stack);
		
		selectWindow("cyst"+toString(nc)+".tif");
		run("Make Substack...", "slices="+zframe+" frames="+i); // Find and duplicate the focused slide
	}
	
	selectWindow("cyst"+toString(nc)+".tif");
	close();
	run("Concatenate...", "all_open"); // Concatenate focused slide in the correct chronological order
	saveAs("Tiff", directory+File.separator+"focused_cysts"+File.separator+"cyst"+toString(nc)+"_focused.tif");
	close();
	close("cyst"+cystnum+".csv");
	
}

// use the focused stack from the previous step to binarize the images for area calculation
File.makeDirectory(directory+File.separator+"binarized_cysts");

for (nc=1; nc<=numCysts; nc++) {

	cystnum = toString(nc);
	
	open(toString(directory+"/focused_cysts/cyst"+toString(nc)+"_focused.tif")); 

	run("Enhance Contrast...", "saturated=0.5 normalize process_all");
	run("Subtract Background...", "rolling=1 create disable stack");
	run("Find Edges", "stack");
	run("Gaussian Blur...", "sigma=2 scaled stack"); // Blur image for binarization
	setAutoThreshold("Default dark");
	
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Default background=Dark");
	
	run("Dilate", "stack"); // Dilate binarized image so the edges close
	run("Dilate", "stack");
	run("Dilate", "stack");
		run("Fill Holes", "stack");
	run("Erode", "stack") ;
	run("Erode", "stack") ;
	run("Erode", "stack") ; // Erode binarized image to counteract the Dilate step

	run("Analyze Particles...", "size=1000-Infinity show=Masks display clear stack");
	
	selectWindow("Mask of cyst"+toString(nc)+"_focused.tif");
	saveAs("Tiff", directory+File.separator+"binarized_cysts"+File.separator+"cyst"+toString(nc)+"_binarized.tif");
	close();
}

  macro "Close All Windows" { 
      while (nImages>0) {
          selectImage(nImages);
          close(); 
      } 
  } 