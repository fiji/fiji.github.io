/*
 * Ths macro will prompt for an OIB file and will segment the cells and run the 3D particle 
 * analyzer on them.
 * 
 * inputs: an appropriate OIB file
 * 		   the corresponding dark image file
 * 
 * output: segmented medium_channel
 *         dark image file (averaged)
 *         estimated background nonuniformity image
 *         statistics table of 3D analyzer
 *         
 * Usage: The user must select an input file (OIB). The user will be prompted to select the
 *        lowest valid slice, and a slice fom whc a background nonuniformity estimate can be made.
 *         
 *    
 *  V3: add background nonuniformity correction       
 *  
 *  V4: changed the thresholder to use a 2-threshold max entropy thresholder. 
 *      the high and low thresholds from the multi-thresholder are used in a
 *      hysteresis thresholder to threshold the cells.
 *      
 *  V5: The threshold is determined as the mean gray level of the object in a pleane that is 2 slices below the center slice time 3x the stdDev of those objects.
 *      This threshold is then  applied to all of the slices.
 *         
 *                   
 * Author: Aryeh Weiss
 * Last modified: 24 June 2014
 * selectImage(npTitle);
 */

VERSION="V5";



MEDIUM_CHANNEL = " - C=0";
NP_CHANNEL = " - C=1";
TX_CHANNEL = " - C=2";

DEFAULT_LOWER_SLICE_DELTA = 2;
DEFAULT_UPPER_SLICE_DELTA = 0;

run("Set Measurements...", "area mean standard modal min center perimeter bounding fit shape stack display redirect=None decimal=3");

// close image windows
run("Close All");	

// test if Log is open. If so, close it.
if (isOpen("Log")) {
	selectWindow("Log");
	run("Close");
}

// test if Summary is open. If so, close it.
if (isOpen("Summary")) {
	selectWindow("Summary");
	run("Close");
}

// test if  Results is open. If so, close it.
if (isOpen("Results")) {
	selectWindow("Results");
	run("Close");
}

// test if  ROI Manager is open. If so, close it.
if (isOpen("ROI Manager")) {
	selectWindow("ROI Manager");
	run("Close");
}

// prompt for input image file. OIB files are assumed.
inputPath = File.openDialog("choose an image file");
inputName = File.getName(inputPath);
inputDir = File.getParent(inputPath);
print("inputPath", inputPath);
print(inputDir);
print(inputName);

/*
fileList = getFileList(inputDirselectImage(npTitle););

for (i=0; i < fileList.length; i++) {
	if (indexOf(fileList[i], "well2") > 0) print(fileList[i]); 
}
*/

inputPrefix = substring(inputName, 0, lastIndexOf(inputName, "."));
print(inputPrefix);

// here we create a uniquely named output directory
// the name includes the inputPrefix, _dapiRoiSetV4, and a date/time string 
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
outputDir = inputPrefix + "_"+VERSION+"_output_"+d2s(year,0)+IJ.pad(d2s(month+1,0),2)+IJ.pad(d2s(dayOfMonth,0),2)+IJ.pad(d2s(hour,0),2)+IJ.pad(d2s(minute,0),2);
outputPath = inputDir + File.separator + outputDir;

// create the output directory. It should never already exists,
// but it is good practice to check
if (!File.exists(outputPath)) File.makeDirectory(outputPath); 

// open dark image
  darkName = "not found";
  list = getFileList(inputDir);
    for (i=0; i<list.length; i++) {
        if (indexOf(list[i], "darkImage") >= 0) {
           	darkName = list[i];
           	darkPrefix = substring(darkName, 0, lastIndexOf(inputName, "."));
           	
           	darkPath = inputDir + File.separator + darkName;
          	run("Bio-Formats Importer", "open=&darkPath autoscale color_mode=Default split_channels view=Hyperstack stack_order=XYCZT");
          	darkTitle = getTitle();
          	run("Grays");
          	
          	selectImage(darkName+TX_CHANNEL);
			close();
			selectImage(darkName+MEDIUM_CHANNEL);
			close();

			selectImage(darkName + NP_CHANNEL);
			rename(inputPrefix+"_NPdark");
			npDarkTitle=getTitle();
			run("Z Project...", "projection=[Average Intensity]");
			npAvgDarkTitle = getTitle;
			run("Median...", "radius=5");
			close(npDarkTitle);

        }
     }

    if (matches(darkName,"not found" )){
  		exit("error: no dark image found!");
    } 
    
// open the image to whic23545 h the *variable* inputPath points.
// rename the three channels to intuitive names
run("Bio-Formats Importer", "open=&inputPath autoscale color_mode=Default split_channels view=Hyperstack stack_order=XYCZT");

// hack to solve bad calibration data written to some OIB files.
getDimensions(width, height, channels, slices, frames);

selectImage(inputName+TX_CHANNEL);
rename(inputPrefix+"_DIC");
dicTitle=getTitle();
close();

selectImage(inputName + NP_CHANNEL);
rename(inputPrefix+"_NP");
npTitle=getTitle();
run("Grays");
run("Properties...", "channels=&channels slices=&slices frames=&frames unit=micron pixel_width=0.207 pixel_height=0.207 voxel_depth=2.0000 frame=[1 sec]");

//close();
selectImage(inputName+MEDIUM_CHANNEL);
rename(inputPrefix+"_medium");
mediumTitle=getTitle();
run("Grays");
run("Properties...", "channels=&channels slices=&slices frames=&frames unit=micron pixel_width=0.207 pixel_height=0.207 voxel_depth=2.0000 frame=[1 sec]");


selectImage(npTitle);
run("Median...", "radius=1 stack");
imageCalculator("Subtract stack", npTitle, npAvgDarkTitle);

selectImage(mediumTitle);
run("Median...", "radius=1 stack");
imageCalculator("Subtract stack", mediumTitle,npAvgDarkTitle);


// insert function to correct background
bgNormTitle = findBgImageEstimate(mediumTitle);


run("Conversions...", " ");
selectImage(mediumTitle);
run("32-bit");
imageCalculator("Divide 32-bit stack", mediumTitle,bgNormTitle);
run("16-bit");



selectImage(npTitle);
run("32-bit");
imageCalculator("Divide 32-bit stack", npTitle,bgNormTitle);
run("16-bit");


//find a slice in the middle which will provide a good threshold in the autothresholder

selectImage(npTitle);
run("Enhance Contrast", "saturated=0.05");
waitForUser("select the lowest valid slice");
selectImage(npTitle);
lowerSlice = getSliceNumber();
print("First slice = ", lowerSlice);


if (lowerSlice > 1) { 
	lowerSlice = lowerSlice - 1;
	selectImage(npTitle);
	run("Slice Remover", "first=1 last=&lowerSlice increment=1");

	selectImage(mediumTitle);
	run("Slice Remover", "first=1 last=&lowerSlice increment=1");

}

getDimensions(width, height, channels, slices, frames);

selectImage(npTitle);
saveAs("zip", outputPath + File.separator + npTitle + ".zip");
rename(npTitle);
npTitle = getTitle();


selectImage(mediumTitle);
saveAs("zip", outputPath + File.separator + mediumTitle + ".zip");
rename(mediumTitle);
mediumTitle = getTitle();

/*
// insert code to prompt for start and end slices 
//getDimensions(width, height, channels, slices, frames);
// prompt for some processing parameters
 Dialog.create("Selection of substack limits for processing");
 
  Dialog.addNumber("Lower slice", DEFAULT_LOWER_SLICE_DELTA);
  Dialog.addNumber("Upper slice", slices - DEFAULT_UPPER_SLICE_DELTA);

//  Dialog.addCheckbox("Ramp", true);
  Dialog.show();
  lowerSlice = Dialog.getNumber();
  upperSlice = Dialog.getNumber();

print("Lower slice: ", lowerSlice);
print("Upper slice: ", upperSlice);

*/

/*
run("Duplicate...", "title="+mediumTitle+"_flatFieldCorrected duplicate");
mediumFlatFieldTitle = getTitle();

blurParam = width/5;
run("Pseudo flat field correction", "blurring=&blurParam stack");
selectImage(mediumFlatFieldTitle+"_background");
backgroundTitle = getTitle();
selectImage(mediumFlatFieldTitle);

*/
// print(width, height, channels, slices, frames);
//find a slice in the middle which will provide a good threshold in the autothresholder

selectImage(mediumTitle);
setSlice(slices/2 - 2);  // middle slice - 2
run("Duplicate...", "title="+mediumTitle+"_bgSlice");
bgSliceTitle = getTitle();

threshold = findThreshold(bgSliceTitle);
selectImage(mediumTitle);
setThreshold(0,threshold);
run("Convert to Mask", "background=Light black");


//setAutoThreshold("Otsu");
/*
thresholds = findMaxEntropyThresholds(bgSliceTitle);
lowThresh = thresholds[0];
highThresh = thresholds[1];
print("low threshold: ",lowThresh, "high threshold: ",highThresh);

close(bgSliceTitle);

selectImage(mediumTitle);
run("Conversions...", "scale");
run("8-bit");
run("3D Hysteresis Thresholding", "high="+d2s(highThresh,0)+" low="+d2s(lowThresh,0));



run("Invert", "stack");


*/
//run("Convert to Mask", "method=Otsu background=Light black");

run("BinaryFilterReconstruct ", "erosions=5 white");
run("Fill Holes", "stack");
run("Watershed", "stack");

run("Analyze Particles...", "size=1000-Infinity pixel show=Masks display exclude clear include summarize add stack");
rename(inputPrefix+"_edgeObjectsRemoved");
maskTitle = getTitle();
run("Grays");

saveAs("zip", outputPath + File.separator + inputPrefix + "_cells.zip");
maskTitle = getTitle();

selectImage(npAvgDarkTitle);
saveAs("zip", outputPath + File.separator + darkPrefix + "_NP.zip");
npAvgDarkTitle = getTitle();
close(npAvgDarkTitle);

selectImage(bgNormTitle);
saveAs("zip", outputPath + File.separator + bgNormTitle + ".zip");
bgNormTitle = getTitle();
close(bgNormTitle);



//close(mediumTitle);
//close(mediumFlatFieldTitle);

selectImage(npTitle);
rollingBallRadius = width/5;
run("Subtract Background...", "rolling=&rollingBallRadius stack");
run("Enhance Contrast", "saturated=0.05");

run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface bounding_box dots_size=5 font_size=10 show_numbers white_numbers store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=&npTitle");

selectImage(maskTitle);
run("3D Objects Counter", "threshold=1 slice=1 min.=1000 max.=3670016 objects surfaces statistics summary");
surfaceTitle = "Surface map of "+ maskTitle + " redirect to " + npTitle;
	initTime = getTime();
	oldTime = initTime;
	while (!isOpen(surfaceTitle)) {
		elapsedTime = getTime() - initTime;
		newTime = getTime() - oldTime;
		if (newTime > 10000) {
			exit("Waited over 10 seconds for 3D OC !");
/*			oldTime = getTime();
			newTime = 0;
			print(elapsedTime/1000, " seconds elapsed"); 
*/			
		}
	
	}
wait(1000);

selectImage("Surface map of "+ maskTitle + " redirect to " + npTitle);
close();

selectImage("Objects map of "+ maskTitle + " redirect to " + npTitle);
saveAs("zip", outputPath + File.separator + inputPrefix + "_3DobjectsMap.zip");

selectWindow("Statistics for "+maskTitle + " redirect to " + npTitle);
saveAs("text", outputPath + File.separator +inputPrefix+"_stats.xls");
close("Results");

setBatchMode(false);

exit("normal exit");


function findBgImageEstimate(inputImageTitle)    {
	//find a slice in the middle which will provide a good threshold in the autothresholder
	waitForUser("select a slice for background estimation");
	bgEstSlice = getSliceNumber();
	print("bgEstSlice = ", bgEstSlice);

	run("Slice Keeper", "first=&bgEstSlice last=&bgEstSlice increment=1");
	rename("bgEstSlice");
	bgEstSliceTitle = getTitle();

	getDimensions(width, height, channels, slices, frames);
	run("Median...", "radius=2");
	
	run("Duplicate...", "title="+bgEstSliceTitle+"_dup");
	bgEstSliceDupTitle = getTitle();
	
	run("Pseudo flat field correction", "blurring=100");
	close(bgEstSliceDupTitle+"_background");
	selectImage(bgEstSliceDupTitle);
	
	setAutoThreshold("Default dark");
	setOption("BlackBackground", true);
	run("Convert to Mask","background=Dark black");
	gridArea=width*height/500;
	run("Grid ", "grid=Lines area=&gridArea color=Black");
	run("Flatten");
	dicedImageTitle = getTitle();
	run("8-bit");
	setThreshold(255, 255);
	run("Convert to Mask", "background=Dark black");
	run("Options...", "iterations=5 count=1 black edm=Overwrite do=Erode");
	run("Analyze Particles...", "  add");
	selectImage(bgEstSliceTitle);
	run("Nonuniform Background Removal", "surface=[Fit Cubic Surface] show");
	close(bgEstSliceTitle+" - Corrected");
	
	rename(mediumTitle+"_bgEst");
	bgNormTitle = getTitle();
	getStatistics(area, mean, min, max, std, histogram);
	run("Conversions...", " ");
	run("32-bit");
	run("Divide...", "value=&max"); // normalize bg image
	
	close(dicedImageTitle);
	close(bgEstSliceTitle);
	close(bgEstSliceDupTitle);

	return bgNormTitle;
}	

function findThreshold(sliceTitle) {
	selectImage(sliceTitle);
	run("Duplicate...", "title="+sliceTitle+"_dup");
	dupTitle = getTitle();

	setAutoThreshold("Otsu");

	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Options...", "iterations=5 count=1 black edm=Overwrite do=Erode");

	run("Set Measurements...", "area mean standard modal min center perimeter bounding fit shape stack display redirect=&sliceTitle decimal=3");

	run("Analyze Particles...", "size=100-40000 pixel circularity=0.0-1.00 display exclude clear include summarize add");

	meanTmp = newArray(nResults);
	stdTmp = newArray(nResults);
	for (i = 0; i < nResults; i++) {
		meanTmp[i] = getResult("Mean", i);
		stdTmp[i] = getResult("StdDev", i);
	}

 	Array.getStatistics(meanTmp, min, max, mean, stdDev);
 	Array.getStatistics(stdTmp, minStd, maxStd, meanStd, stdDev);
 	close(dupTitle);
 
	return mean+(3*meanStd);
}


function findMaxEntropyThresholds(inputTitle) {

/* 
 *  This function find two threhsolds for an input image, using 
 *  a max entropy mutli-thresholding plugin. The thresholds are returned 
 *  in an array called thresholds, arranged as (lowThresh, highThresh)
 *
 */
	run("Duplicate...", "title="+inputTitle+"_dup");
	inputDupTitle = getTitle();

	
	run("Conversions...", "scale");
	run("8-bit");
	run("Maximum Entropy Multi-Threshold", "number=2");

	selectWindow("Log");
	logStringArray = split(getInfo("window.contents"),'\n');
	tmpString = logStringArray[logStringArray.length-2];
//	print(logStringArray[logStringArray.length-2]);
	tmpStringArray = split(tmpString, " ");

	lowThresh = d2s(tmpStringArray[tmpStringArray.length-2],0);
	highThresh = d2s(tmpStringArray[tmpStringArray.length-1],0);
	print("low threshold: ",lowThresh, "high threshold: ",highThresh);

	close(inputDupTitle);
	thresholds = newArray(lowThresh, highThresh);
	return thresholds;
}
