

//The purpose of this macro is to analyze a folders worth of in vitro condensate data.
//this macro outputs several pieces of data. it outputs the size and the shape of each condensate. 
//it also outputs the brightness of a center of each ROI in both channels.
//we use the centroid to avoid an interface effect and more accuratly measure the dense phase

//from each input file you will get 3 output files representing the RNA intensity, the protein intensity and the 
//shape and size parameters. 

//asks for directory
input = getDir("Choose input directory")
output = getDir("Choose output directory")




//we don't use any really tiny condensates because they seem more likely 
//to be junk or hot pixels 
//this can be adjusted to a different minimum size of a blob here
pixel_size = 9


//setBatchMode(true);
list = getFileList(input);
for (i = 0; i < list.length; i++)
    getBlobStats(input, output, list[i]);
   
   print("Done");
//this just adds some words to the end of the log so that I can make sure it didn't glich out
//part way through the run
	
   



function getBlobStats(input, output, filename) {

	//open current file
	run("Bio-Formats", "open=" + input + filename);

	
	//define variables for slice logic
	currentarea = 0;
	maxarea = 1;
	myslice = 0;


	//Determines wavelengths from metadata
	infoString = getMetadata("Info");
	titleString = getTitle();
	getDimensions(width, height, channels, slices, frames);



	
	name = File.nameWithoutExtension;
	progood = 1; 


	
	ProChannel = 1;  //in my images protein is the first channel adjust as necessary
	RNAChan = 2 ;    //change as necessary
					// all of the thresholding is done on a single channel and then used for both channels
					//in my images the protein channel has higher SNR but this could be different

    if (maxarea > 0 && progood == 1) {
    	
		Stack.setChannel(ProChannel);
		string=getInfo("log");


		//this block duplicates the original into a number of smaller files
		//one is for setting the ROIs and the other are to measure each channel
			
		selectWindow(filename);
		run("Duplicate...", "title=onepro duplicate channels=ProChannel");  //duplicates the protein channel 

		
		selectWindow("onepro");
		run("Z Project...", "projection=[Max Intensity]"); //z projects the channel
		
		
		selectWindow("MAX_onepro");
		run("Duplicate...", "title=MAX_two duplicate");

		selectWindow("onepro"); // closes the multislice so we only have the z-projection
		run("Close");
		selectWindow("MAX_onepro");

		Stack.setChannel(RNAChan);
		selectWindow(filename);
		run("Duplicate...", "title=onerna duplicate channels=RNAChan");

		selectWindow("onerna");
		run("Z Project...", "projection=[Max Intensity]"); //z projects the channel
		selectWindow("onerna"); // closes the multislice so we only have the z-projection
		run("Close");

		selectWindow("MAX_two");
		
		//sets the type of threshold to use on all the images

			
		setAutoThreshold("Yen dark");  //type of threshold seems to matter. its probably worth trying other options
		setOption("BlackBackground", true);
		run("Convert to Mask");  //applies the threshold
		
		run("Watershed");  // runs a watershed 
							// we use a watershed to segment connected droplets into individuals for
							//measurment purposes. 
							//this could be adjusted

		
		

		run("Set Measurements...", "area centroid shape decimal=3"); //first measurments are just shape parameters on
																		// the full image


		run("Analyze Particles...", "size=9-Infinity pixel exclude clear include add");  //analyzes particles and puts them in ROI manager
		
		selectWindow("ROI Manager");
		roiManager("Measure");


		if(nResults() > 0){ //saves the results if any condensates were measured. 
		
		selectWindow("Results");
		saveAs("Results", output + name + "_area.tsv");
		

		if (isOpen("Results")) { 
			selectWindow("Results"); 
         	run("Close");
    	} 

    
    
    } else {
    	print("no area to measure");
    	string=getInfo("log");
		//File.saveString(string, output + date + "_log.txt");
    	
    }
		
		
			//we wanted to measure only the centroid of the condensates based on some work we did 
			//that found a substantial interface
			//this bit of code is very kludgy (and slow) but it takes the existing ROIs and one at a time converts
			//them to a 3,3 box centered on the original ROI. 
			//this is used for the intensity measurments.	
		
		selectWindow("ROI Manager");
		count = roiManager("count");  //figures out how many particles there are.
		
		for (i=0; i<count; i++){   
		 	roiManager("Select", 0);
			run("Specify...", "width=3 height=3 centered");
			roiManager("Add");
			roiManager("select", 0);
			roiManager("delete");
		 }

		 
		selectWindow("MAX_onepro");
		
		run("Set Measurements...", "mean redirect=None decimal=3"); //only measures the intensity
		 
		roiManager("Show All");
		
		roiManager("Measure");

		
	
	if(nResults() > 0){
		
		selectWindow("Results");

			
		saveAs("Results", output + name + "_protein.tsv");
		

		if (isOpen("Results")) { 
			selectWindow("Results"); 
         	run("Close");
    	} 

    
    
    } else {
    	print("no area to measure");
    	string=getInfo("log");
		File.saveString(string, output +  "_log.txt");
    	
    }

    
		selectWindow("MAX_onerna"); //now we do the same measures on the RNA channel.
		
		run("Set Measurements...", "mean redirect=None decimal=3");
		 
		roiManager("Show All");
		
		roiManager("Measure");

		
	
	if(nResults() > 0){
		
		selectWindow("Results");

			
		saveAs("Results", output + name + "_rna.tsv");
		

		if (isOpen("Results")) { 
			selectWindow("Results"); 
         	run("Close");
    	} 

    
    
    } else {
    	print("no area to measure");
    	string=getInfo("log");
		File.saveString(string, output +  "_log.txt");
    	
    }
	
    // close duplicate window
    	selectWindow("MAX_onepro");
    	run("Close");
    	selectWindow("MAX_two");
    	run("Close");
    	selectWindow("MAX_onerna");
    	run("Close");
    	roiManager("reset");
    	close("*");
    	
    }
}


// function returning the number of the column which name is contained in kind. return -1 if doesn't exists
function columnLabelList (kind) {

	columnNumber=-1;
	if (nResults() > 0 && isOpen("Results")) {
		selectWindow("Results");
   		results = getInfo();
   		lines = split(results, "\n");
  		headings = lines[0];
		titlesofcolumns = split(headings, ",\t");
		for (a=0; a<titlesofcolumns.length; a++) {if (titlesofcolumns[a] == kind) columnNumber=a;}
	}
	return columnNumber;
}



// function getting the sum of the value from the column name contained in kind
function sumColumnTab (kind) {
	sum = 0;
	if (columnLabelList (kind) >=0) {
		for (a=0; a<nResults(); a++) {
			sum = sum + getResult(kind,a);
		}
	return sum;
}

