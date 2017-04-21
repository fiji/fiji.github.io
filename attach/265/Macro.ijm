dir = getDirectory("Choose Sorce Directory"); //with slash at the end
dir2 = getDirectory("Choose Destination Directory");
acq = getFileList(dir); //this is an array of dir starting from 0 and with slash at the end
setBatchMode(true); 
for (i=0; i<acq.length; i++) {
	if (endsWith(acq[i], "tif")) {
		open(dir+acq[i]);
		run("Split Channels");
		close("C3-"+acq[i]);
		imageCalculator("Subtract stack", "C2-"+acq[i], "C1-"+acq[i]);
		//Save Z-projection of 169
		selectWindow("C2-"+acq[i]);
		run("Despeckle");
		run("Z Project...", "projection=[Average Intensity]");
		selectWindow("AVG_C2-"+acq[i]);
		run("Save", "save=["+dir2+replace(acq[i],"/","-Zproj-169-BGsubtr.tif]"));
		close("AVG_C2-"+acq[i]);
		//Save Z-projection of SHG
		selectWindow("C4-"+acq[i]);
		run("Despeckle");
		run("Z Project...", "projection=[Average Intensity]");
		selectWindow("AVG_C4-"+acq[i]);
		run("Save", "save=["+dir2+replace(acq[i],"/","-Zproj-SHG-BGsubtr.tif]"));
		close("*");
	}
}
	