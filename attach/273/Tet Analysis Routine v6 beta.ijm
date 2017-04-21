setBatchMode(false);
 /*
 * This macro code will take a multi-channel TIFF stack of a Tetrahymena cell stained with a robust basal body marker
 * (i.e., Centrin) and identify the location of all basal bodies within the image as well as each basal body's most 
 * logical partner within a ciliary row. The code is split up into 
 * 
 */
initMacroReturn=initMacro();//Initializes the macro by getting the stage of the macro that needs to be run and the folder structure for what the macro should process.
/*
 * Part 1-Find BBs
 */
if(initMacroReturn[0]=="1.  Find BBs"){ 
	initPart1Return=initPart1(initMacroReturn[1]);//Initializes the macro by getting the: 1) channel names 2) the channel number for the robust BB marker of interest and 3) whether the macro is to be run on a single image or a directory of images. This information is returned as an array. 
	fileList=createFileList(initMacroReturn[1],initPart1Return[7]);//Creates a list of all the files (directory+file name) that the macro will run on. Returns the file list as an array.
	for(i=0; i<fileList.length; i++){//The start of the loop that marches through the file list.
		saveDir=File.getParent(fileList[i])+File.separator;//Assigns the parent directory of the image file to the variable saveDir.
		saveName=File.getName(fileList[i]);//Assigns the file name of the image file to the variable saveName.
		saveNameNoExt=substring(saveName,0,lengthOf(saveName)-4);//Removes the last 4 characters of the file name (i.e., ".tif") and assigns the extensionless name to the variable saveNameNoExt.
		writeTo=File.open(saveDir+saveNameNoExt+".txt");//Creates a .txt file in the parent directory with the same name as the image file. All image information will be written to this file. The file is assigned the varaible writeTo.
		writeToOpenFile(writeTo,initPart1Return,"Dialog");//Writes the dialog information from line 5 to the writeTo file. writeToOpenFile is a custom function that adds text to an open .txt file.
		title=fileList[i];//The bio-formats importer will not accept an array entry as an input, so the file (directory+file name) is assigned to the variable title.
		run("Bio-Formats Importer", "open=&title autoscale color_mode=Default split_channels view=Hyperstack stack_order=XYCZT");//Opens the file using bio formats which opens .tifs. nd2s and other propriety files. 
		getVoxelSize(voxWidth,voxHeight,voxDepth,unit);//Returns the voxel size of the opened image. 
		if(voxWidth<0.1){//Start conditional IF statement for images with voxels less than 0.1 microns
			for(ii=0; ii<nImages; ii++){//Start loop reducing image size
				selectImage(ii+1);//Select the image.
				width=getWidth;//Get the image width.
				height=getHeight;//Get the image height.
				depth=nSlices;//Get the image depth.
				newWidth=round(width/2);//Divides the width by 2 and returns the value as the variable newWidth.
				newHeight=round(height/2);//Divides the height by 2 and returns the value as the variable newHeight.
				run("Size...", "width=&newWidth height=&newHeight depth=&depth average interpolation=Bicubic");//The imageJ function for reducing the size of the image.
			}//End loop reducing image size
		}//End conditional IF statement for images with voxels less than 0.1 microns
		selectImage(parseFloat(initPart1Return[1]));//dialog[1] is the value that the user selected as the image that contains the robust BB marker. Therefore, this activates the image with the robust BB marker.
		features=roughFeatureExtraction(1, 1);//roughFeatureExraction is a function that scans a max intensity projection for large objects that appear after large radius blurring. The function returns the XY coordinates for outlines encompassing each feature.
		shapeFilteredFeatures=shapeFilterRois(features,1,1);//shapeFilteredFeatures is a function that deletes features that are unlikely to be Tetrahymena from the list of features generated from roughFeatureExtraction. The function returns the XY coordinates of the most Tetrahymena-like features.
		largestFeature=returnLargestRoi(shapeFilteredFeatures,200);//returnLargestRoi is a function that scans a list of features and simply returns the XY coordinates of the largest feature. This feature is expected to be the largest Tetrahymena cell in the image stack
		writeToOpenFile(writeTo,largestFeature,"Largest Feature");//Write the coordinates of the largest Tetrahymena cell to the writeTo file.
		croppedFeature=cropXy(largestFeature,5,1);//croppedFeature is a function that crops an image using a feature of interest. It returns the coordinates of the feature within the resulting cropped image, which is necessary because the absolute coordinates for each part of the selection will be different after cropping.
		writeToOpenFile(writeTo,croppedFeature,"Cropped Feature");//Write the coordinates of the cropped largest Tetrahymena cell to the writeTo file.
		selectImage(parseFloat(initPart1Return[1]));//dialog[1] is the value that the user selected as the image that contains the robust BB marker. Therefore, this activates the image with the robust BB marker.
		coords=findLocalMaxima(croppedFeature,6);//findLocalMaxima is a function that returns the XYZ coordinates of the maximum pixel intensity of each BB in the image. The function returns a two part array (3(XYZ)/1(intensity)) with the last part being intensity of each object and the first part being the xyz coords of each object
		points=cleanUpCoordArray(coords);//cleanUpCoordArray is a function that removes the strings "X", "Y", and "Z" from the array generated by findLocalMaxima.
		convertArrayToRois(coords,1);//convertArrayToRois is a function that takes a "clean" array of coordinates (XYZ) and converts them to regions of interest usable by the ImageJ region of interest manager.
		bbsAndOa=skewKurtFilterHomog(4);//skewKurtFilterHomog is a function that examines the local intensity environment of a region of interest. In this instance, it will remove regions of interest (i.e., potential BBs) whose local skew and kurtosis are low or negative which indicates a homogenous intensity background. Since the OA is characterized by such BBs, this function also finds the OA (returned as a XYZ of the center of mass for homogenous objects). It returns the XYZ coordinates of all objects that were not rejected due to skew and kurtosis.
		bbs=Array.slice(bbsAndOa,0,bbsAndOa.length-3);//Slice off the portion of bbsAndOA array that contains the actual BBs.
		oa=Array.slice(bbsAndOa,bbsAndOa.length-3,bbsAndOa.length);//Slice off the portion of bbsAndOA array that contains the center of mass for the OA.
		writeToOpenFile(writeTo,oa,"OA");//Writes the coordinates of the OA center of mass to the writeTo file.
		waitForUser("1");
		removeMidBbs(2.25);//removeMidBbs is a function that examines the convex hull of the true BBs and looks for any true BBs that further than some distance from the surface of the convex hull. The function removes these BBs from the list as they are internalized BBs and should not be used for cellular geometry.
		finalBbs=convertRoisToCoords(1);//convertRoisToCoords is a function that converts regions of interest from the region of interest manager into XYZ coordinates in the form of an array.
		writeToOpenFile(writeTo,finalBbs,"Final BBs");//Writes the coordinates of the final BBs to the writeTo file.
		waitForUser("2");
		poles=findPoles(1,bbs,5);//findPoles is a function that defines the anterior and posterior pole by finding the 5 pairs of BBs that have the greatest separation distance. Since Tet are ellipsoids the location of the pairs that have the greatest separtion distance reflects the cell poles. The anterior pole is defined by exmaining the integrated intensity of the area surrounding the pole, which will have a greater integrated intensity due to the OA and the closely spaced BBs at the tip.
		writeToOpenFile(writeTo,poles,"Poles");//Writes the coordinates of the poles to the writeTo file.
		File.close(writeTo);//The write to file must be closed since only a single writeTo file can be open at any time
		run("Close All");//Closes all the images that are currently open
		roiManager("Reset");//Resets the ROI manager so it is blank before the next file is run
	}
}
/*
 * Part 2-Quality Control BBs
 */
if(initMacroReturn[0]=="2.  Quality Control BBs"){
	initPart2Return=initPart2(initMacroReturn[1]);//Initializes Part 2 of the macro
	fileList=createFileList(initMacroReturn[1],initPart2Return[0]);//Creates a list of all the files (directory+file name) that the macro will run on. Returns the file list as an array.
	for(i=0; i<fileList.length; i++){//The start of the loop that marches through the file list.
		saveDir=File.getParent(fileList[i])+File.separator;//Assigns the parent directory of the image file to the variable saveDir.
		saveName=File.getName(fileList[i]);//Assigns the file name of the image file to the variable saveName.
		saveNameNoExt=substring(saveName,0,lengthOf(saveName)-4);//Removes the last 4 characters of the file name (i.e., ".tif") and assigns the extensionless name to the variable saveNameNoExt.
		textFile=saveDir+saveNameNoExt+".txt";//Assigns the file name of the text file that was created for each image in Part 1 to the variable textFile
		textFileArray=openTextFile(textFile);//Opens the text file and returns its contents as an array named textFileArray
		dialogText=Array.slice(textFileArray,0,9);//Positions 0-10 of the textFileArray contain the dialog information that was gathered from the user in Part 1.
		cellOutline=grabLargestFeature(textFileArray,"Largest Feature","Cropped Feature");
		cropBox=grabCroppedFeature(textFileArray,"Cropped Feature","OA");
		oaFeature=grabOaFeature(textFileArray,"OA","Final BBs");
		bBsFeature=grabBbFeature(textFileArray,"Final BBs","Poles");
		poleFeature=grabPoleFeature(textFileArray,"Poles");
		title=fileList[i];//The bio-formats importer will not accept an array entry as an input, so the file (directory+file name) is assigned to the variable title.
		run("Bio-Formats Importer", "open=&title autoscale color_mode=Default split_channels view=Hyperstack stack_order=XYCZT");//Opens the file using bio formats which opens .tifs. nd2s and other propriety files. 
		getVoxelSize(voxW,voxH,voxD,unit);
		for(ii=0; ii<nImages; ii++){
			selectImage(ii+1);
			if(ii+1==parseFloat(dialogText[2])){
				bbImage=getTitle;
			}
			createSelectionAndCrop(cellOutline,5);//This will take the cellOutline and enlarge it by 5 microns
		}
		convertArrayToRois(bBsFeature,1);
		bbCentroids=convertRoiToBinary();
		run("16-bit");
		run("Multiply...", "value=65535 stack");
		bbCentroidsTitle=getTitle;
		/*
		 * Remove extraneous BBs
		 */
		run("Merge Channels...", "c1=&bbImage c2=&bbCentroidsTitle create keep");
		mergeImage=getImageID();
		setBatchMode("show");
		deleteBBs=interactiveCursor(1,2,mergeImage,bbCentroids);//input for width and depth are in microns
		setBatchMode("hide");
		if(deleteBBs.length>0){
			modifyBinaryImage(bbCentroidsTitle,"delete",deleteBBs,1,"square");
			selectImage(mergeImage);
			close();
		}
		/*
		 * Add missed BBs
		 */
		run("Merge Channels...", "c1=&bbImage c2=&bbCentroidsTitle create keep");
		mergeImage=getImageID();
		setBatchMode("show");
		addBBs=interactiveCursor(1,2,mergeImage,bbImage);//input for width and depth are in microns
		setBatchMode("hide");
		if(addBBs.length>0){
			modifyBinaryImage(bbCentroidsTitle,"add",addBBs,1,"square");
			selectImage(mergeImage);
			close();
		}
		 /*
		  * Rewrite the text file if the user updated the BBs 
		  */
		if(deleteBBs.length>0 || addBBs.length>0){
			roiManager("Reset");
			selectImage(bbCentroidsTitle);
			makeBinaryImage();
			bbs=convertSelectionsToRois(0,1);
			finalBbsCoords=convertRoisToCoords(1);
			File.delete(saveDir+saveNameNoExt+".txt");
			writeTo=File.open(saveDir+saveNameNoExt+".txt");//Creates a .txt file in the parent directory with the same name as the image file. All image information will be written to this file. The file is assigned the varaible writeTo.			
			writeToOpenFile(writeTo,dialogText,"None");//Writes the dialog information from line 5 to the writeTo file. writeToOpenFile is a custom function that adds text to an open .txt file.
			writeToOpenFile(writeTo,cellOutline,"Largest Feature");
			writeToOpenFile(writeTo,cropBox,"Cropped Feature");
			writeToOpenFile(writeTo,oaFeature,"OA");
			writeToOpenFile(writeTo,finalBbsCoords,"Final BBs");
			writeToOpenFile(writeTo,poleFeature,"Poles");//Writes the coordinates of the poles to the writeTo file.
		}
		File.close(writeTo);//The write to file must be closed since only a single writeTo file can be open at any time
		run("Close All");//Closes all the images that are currently open
		roiManager("Reset");//Resets the ROI manager so it is blank before the next file is run
	}
}
/*
 * Part 3-Make Connections
 */
if(initMacroReturn[0]=="3.  Make Connections"){
	initPart2Return=initPart2(initMacroReturn[1]);//Initializes Part 2 of the macro
	fileList=createFileList(initMacroReturn[1],initPart2Return[0]);//Creates a list of all the files (directory+file name) that the macro will run on. Returns the file list as an array.
	for(i=0; i<fileList.length; i++){//The start of the loop that marches through the file list.
		saveDir=File.getParent(fileList[i])+File.separator;//Assigns the parent directory of the image file to the variable saveDir.
		saveName=File.getName(fileList[i]);//Assigns the file name of the image file to the variable saveName.
		saveNameNoExt=substring(saveName,0,lengthOf(saveName)-4);//Removes the last 4 characters of the file name (i.e., ".tif") and assigns the extensionless name to the variable saveNameNoExt.
		textFile=saveDir+saveNameNoExt+".txt";//Assigns the file name of the text file that was created for each image in Part 1 to the variable textFile
		textFileArray=openTextFile(textFile);//Opens the text file and returns its contents as an array named textFileArray
		dialogText=Array.slice(textFileArray,0,9);//Positions 0-10 of the textFileArray contain the dialog information that was gathered from the user in Part 1.
		cellOutline=grabLargestFeature(textFileArray,"Largest Feature","Cropped Feature");
		cropBox=grabCroppedFeature(textFileArray,"Cropped Feature","OA");
		oaFeature=grabOaFeature(textFileArray,"OA","Final BBs");
		bBsFeature=grabBbFeature(textFileArray,"Final BBs","Poles");
		poleFeature=grabPoleFeature(textFileArray,"Poles");
		poleNumbers=parseArray(poleFeature);
		title=fileList[i];//The bio-formats importer will not accept an array entry as an input, so the file (directory+file name) is assigned to the variable title.
		run("Bio-Formats Importer", "open=&title autoscale color_mode=Default split_channels view=Hyperstack stack_order=XYCZT");//Opens the file using bio formats which opens .tifs. nd2s and other propriety files. 
		getVoxelSize(voxW,voxH,voxD,unit);
		for(ii=0; ii<nImages; ii++){
			selectImage(ii+1);
			if(ii+1==parseFloat(dialogText[2])){
				bbImage=getTitle;
			}
			createSelectionAndCrop(cellOutline,5);//This will take the cellOutline and enlarge it by 5 microns
		}	
		convertArrayToRois(bBsFeature,1);
		bbCentroids=convertRoiToBinary();
		run("16-bit");
		run("Multiply...", "value=65535 stack");
		bbCentroidsTitle=getTitle;
		makeBinaryImage();
		roiManager("Reset");
		bbs=convertSelectionsToRois(0,1);
		finalBbsCoords=convertRoisToCoords(1);
		connections=iterativeBasalBodyFinding(bbs,15,poleNumbers);//iterativeBasalBodyFinding is a function that connects each BB to its anterior partner. In a perfect cell without chirality, the anterior partner should be a BB that is coplanar with the BB, the anterior pole and the posterior pole. So this procedure scores BBs based on their deviation from this plane and the distance between the potential anterior BB and the BB itself. The connected BB will minimize the deviation from the plane and the distance from the BB and it must be closer to the anterior pole than the BB itself.
		File.delete(saveDir+saveNameNoExt+".txt");
		writeTo=File.open(saveDir+saveNameNoExt+".txt");//Creates a .txt file in the parent directory with the same name as the image file. All image information will be written to this file. The file is assigned the varaible writeTo.			
		writeToOpenFile(writeTo,dialogText,"None");//Writes the dialog information from line 5 to the writeTo file. writeToOpenFile is a custom function that adds text to an open .txt file.
		writeToOpenFile(writeTo,cellOutline,"Largest Feature");
		writeToOpenFile(writeTo,cropBox,"Cropped Feature");
		writeToOpenFile(writeTo,oaFeature,"OA");
		writeToOpenFile(writeTo,finalBbsCoords,"Final BBs");
		writeToOpenFile(writeTo,poleFeature,"Poles");//Writes the coordinates of the poles to the writeTo file.
		writeToOpenFile(writeTo,connections,"Connections");//Writes the coordinates of the connected BBs to the writeTo file.
		File.close(writeTo);//The write to file must be closed since only a single writeTo file can be open at any time
		run("Close All");//Closes all the images that are currently open
		roiManager("Reset");//Resets the ROI manager so it is blank before the next file is run
}

/*
 * 
 * 
 * 
 * Below are the user-created functions that are utilized by the macro
 * 
 * 
 * 
 * 
 */
/*
 * 
 */
function initMacro(){
	stage=newArray("1.  Find BBs", "2.  Quality Control BBs", "3.  Make Connections", "4.  Quality Control Connections", "5.  Create Visual 3D Models");
	structure=newArray("Single File","Single Folder", "Folder of Folders");//An array to populate the dialog with choices for the type of folder structure to analyze.
 	Dialog.create("Batch Analysis of BB Spatial Distribution");//Creates a dialog called "Batch Analysis of BB Spatial Distribution
 	Dialog.addChoice("Which stage of analysis?", stage);//A dialog line asking the user to define how many channels are in the images in the folder structure to be analyzed
 	Dialog.addChoice("File structure:", structure);
 	Dialog.show();//ImageJ function to actually show the dialog on the screen.
 	stage=Dialog.getChoice;//Assign the dialog input for the total number of channels
 	fileStructure=Dialog.getChoice;
 	initMacroArray=newArray(stage,fileStructure);
	return initMacroArray;
}
/*
 * initPart1 is a function that runs all the dialogs at the start of the macro that gathers information about the analysis.
 */
function initPart1(fileType){//No arguments.
	dialogArray=batchDialog();//The function batchDialog will create a dialog suited for working on an unopened image(s).
	if(fileType=="Single File"){//The 6th item returned by the batchDialog is the file structure chosen by the user. In this case it is a "single file"
		file=File.openDialog("Choose File");//Pops a file window to navigate to image file.
		dialogArray=Array.concat(dialogArray,file);//Adds the path of the file to the dialogArray
	}
	if(fileType=="Single Folder"){//The 6th item returned by the batchDialog is the file structure chosen by the user. In this case it is a "single folder"
		dir=getDirectory("Choose Folder");//Pops a file window to navigate to folder.
		dialogArray=Array.concat(dialogArray,dir);//Adds the path the folder to the dialogArray
	}
	if(fileType=="Folder of Folders"){//The 6th item returned by the batchDialog is the file structure chosen by the user. In this case it is a "folder of folders"
		parent=getDirectory("Choose Parent Folder");//Pops a file window to navigate to folder of folders.
		dialogArray=Array.concat(dialogArray,parent);//Adds the path the folder of folders to the dialogArray
	}
	return dialogArray;//Returns the dialogArray that was just created
}
/*
 * batchDialog is a function that runs if there is an image open when the code is started
 */
 function batchDialog(){//No arguments 	
 	channels=newArray(1,2,3,4);//An array to populate the dialog with choices for which channel contains the robust BB marker.
 	Dialog.create("Batch Analysis of BB Spatial Distribution");//Creates a dialog called "Batch Analysis of BB Spatial Distribution
 	Dialog.addChoice("How many channels are your images?", channels);//A dialog line asking the user to define how many channels are in the images in the folder structure to be analyzed
 	Dialog.addChoice("Which channel is the robust BB marker?", channels);//A dialog line asking the user to pick which channel has the robust BB marker (i.e., which channel should be used for automated BB detection).
 	Dialog.addString("Channel 1", "");//A dialog line waiting for the user to input the name of Channel 1
 	Dialog.addString("Channel 2", "");//A dialog line waiting for the user to input the name of Channel 2
 	Dialog.addString("Channel 3", "");//A dialog line waiting for the user to input the name of Channel 3
 	Dialog.addString("Channel 4", "");//A dialog line waiting for the user to input the name of Channel 4
 	Dialog.show();//ImageJ function to actually show the dialog on the screen.
 	totalChannels=Dialog.getChoice;//Assign the dialog input for the total number of channels
 	robustChannel=Dialog.getChoice;//Assign the dialog input for the channel that contains the robust BB marker
 	channel1=Dialog.getString;//Retrieve the name for channel 1
 	images=0;//The images variable will be incremented to match the number of named channels in the image.
 	if(lengthOf(channel1)>0){//If the user inputed a string into the Channel 1 dialog, the increase images by 1.
 		images++;
 	}
 	channel2=Dialog.getString;//Retrieve the name for channel 2
  	if(lengthOf(channel2)>0){//If the user inputed a string into the Channel 1 dialog, the increase images by 2.
 		images++;
 	}
 	channel3=Dialog.getString;//Retrieve the name for channel 3
  	if(lengthOf(channel3)>0){//If the user inputed a string into the Channel 1 dialog, the increase images by 3.
 		images++;
 	}
 	channel4=Dialog.getString;//Retrieve the name for channel 4
  	if(lengthOf(channel4)>0){//If the user inputed a string into the Channel 1 dialog, the increase images by 4.
 		images++;
 	}
 	dialogArray=newArray(totalChannels,robustChannel,channel1,channel2,channel3,channel4,images);//The dialogArray will be filled with the information gathered from the dialog.
 	return dialogArray;
}
/*
 * Creates a file list for the main part of the macro to process if it is in batch mode.
 */
function createFileList(fileStructure,file){//Arguments 1)"fileStructure" is the file structure that the analysis will work on, 2)"file" is the path to the image, folder or folder of folders, 3) "images" is the number channels for each image in the folder structure
	filesToProcess=newArray();
	if(fileStructure=="Single File"){
		dir=File.getParent(file)+File.separator;
		name=File.getName(file);
		path=dir+name;
		filesToProcess=Array.concat(filesToProcess,path);	
	}
	if(fileStructure=="Single Folder"){
		dir=file;
		temp=getFileList(dir);
		filesToProcess=newArray();
		for(i=0; i<temp.length; i++){
			if (endsWith(temp[i],".nd2")==true || endsWith(temp[i],".tif")==true){
				add=dir+temp[i];
				filesToProcess=Array.concat(filesToProcess,add);
			}
		}
	}
	if(fileStructure=="Folder of Folders"){
		dir=file;
		foldersToProcess=getFileList(dir);
		filesToProcess=newArray();
		for(i=0; i<foldersToProcess.length; i++){
			temp=getFileList(dir+foldersToProcess[i]);
			for(ii=0; ii<temp.length; ii++){
				if (endsWith(temp[ii],".nd2")==true /*|| endsWith(temp[ii],".tif")==true*/){
					add=dir+foldersToProcess[i]+temp[ii];
					filesToProcess=Array.concat(filesToProcess,add);
				}
			}
		}	
	}
	return filesToProcess;
}
/*
 * initPart2 is a function that runs all the dialogs at the start of the macro that gathers information about the analysis.
 */
function initPart2(fileType){//No arguments.
	dialogArray=newArray();
	if(fileType=="Single File"){//The 6th item returned by the batchDialog is the file structure chosen by the user. In this case it is a "single file"
		file=File.openDialog("Choose File");//Pops a file window to navigate to image file.
		dialogArray=Array.concat(dialogArray,file,0);//Adds the path of the file to the dialogArray
	}
	if(fileType=="Single Folder"){//The 6th item returned by the batchDialog is the file structure chosen by the user. In this case it is a "single folder"
		dir=getDirectory("Choose Folder");//Pops a file window to navigate to folder.
		dialogArray=Array.concat(dialogArray,dir,1);//Adds the path the folder to the dialogArray
	}
	if(fileType=="Folder of Folders"){//The 6th item returned by the batchDialog is the file structure chosen by the user. In this case it is a "folder of folders"
		parent=getDirectory("Choose Parent Folder");//Pops a file window to navigate to folder of folders.
		dialogArray=Array.concat(dialogArray,parent,1);//Adds the path the folder of folders to the dialogArray
	}
	return dialogArray;//Returns the dialogArray that was just created
}
/*
 * writeToOpenFile writes to a text file that is generated for each raw image data file. The writeTo file is generated by the File.open() ImageJ macro command. The eventual contents of the writeTo file will be strings not numbers so parseFloat will be needed to convert string numbers into floating point values.
 */
function writeToOpenFile(file,array,feature){//Arguments: 1) "file" is the name of the text file generated by File.open(), 2) "array" is an array of strings or values or both that will be written to the writeTo file with each entry being listed on a new line, 3) "feature" is a string that is the title for the entry into the writeTo file that actually gets written before any of the actual data entries. This is useful because it demarcates where the list of features starts in the writeTo file.
	if(feature!="None"){
		print(file, feature);//Print "feature" to the writeTo file.
	}
	for(i=0; i<array.length; i++){//Start loop for writing each individual element of "array" to the writeTo file.
		print(file,array[i]);//Write the individual element of "array" to the writeTo file.
	}//End loop for writing each individual element of "array" to the writeTo file.
}
/*
 * 
 */
function openTextFile(file){
	parameters=File.openAsString(file);
	parameterArray=split(parameters,"\n");
	return parameterArray;
}
/*
 * 
 */
function grabLargestFeature(parameterArray,heading1,heading2){
	chunkArray=newArray(2);
	for(ii=0; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading1)==true){
			chunkArray[0]=ii+1;
			ii=parameterArray.length;
		}
	}
	for(ii=chunkArray[0]; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading2)==true){
			chunkArray[1]=ii;
			ii=parameterArray.length;
		}
	}
	largestFeature=Array.slice(parameterArray,chunkArray[0],chunkArray[1]);
	return largestFeature;
}
/*
 * 
 */
function grabCroppedFeature(parameterArray,heading1,heading2){
	chunkArray=newArray(2);
	for(ii=0; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading1)==true){
			chunkArray[0]=ii+1;
			ii=parameterArray.length;
		}
	}
	for(ii=chunkArray[0]; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading2)==true){
			chunkArray[1]=ii;
			ii=parameterArray.length;
		}
	}
	croppedFeature=Array.slice(parameterArray,chunkArray[0],chunkArray[1]);
	return croppedFeature;
}
/*
 * 
 */
function grabOaFeature(parameterArray,heading1,heading2){
	chunkArray=newArray(2);
	for(ii=0; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading1)==true){
			chunkArray[0]=ii+1;
		}
	}
	for(ii=chunkArray[0]; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading2)==true){
			chunkArray[1]=ii;
			ii=parameterArray.length; 
		}
	}
	OAFeature=Array.slice(parameterArray,chunkArray[0],chunkArray[1]);
	return OAFeature;
}
/*
 * 
 */
function grabBbFeature(parameterArray,heading1,heading2){
	chunkArray=newArray(2);
	for(ii=0; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading1)==true){
			chunkArray[0]=ii+1;
			ii=parameterArray.length;
		}
	}
	for(ii=chunkArray[0]; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading2)==true){
			chunkArray[1]=ii;
			ii=parameterArray.length; 
		}
	}
	bbFeature=Array.slice(parameterArray,chunkArray[0],chunkArray[1]);
	return bbFeature;
}
/*
 * 
 */
function grabPoleFeature(parameterArray,heading1){
	chunkArray=newArray(2);
	for(ii=0; ii<parameterArray.length; ii++){
		temp=parameterArray[ii];
		if(endsWith(temp,heading1)==true){
			chunkArray[0]=ii+1;
			ii=parameterArray.length;
		}
	}
	poleFeature=Array.slice(parameterArray,chunkArray[0],parameterArray.length);
	return poleFeature;
}
/*
 * 
 */
function createSelectionAndCrop(array,enlarge){
	chunkArray=newArray(4);
	for(i=0; i<array.length; i++){
		if(array[i]=="x"){
			chunkArray[0]=i+1;
		}
		if(array[i]=="y"){
			chunkArray[1]=i;
			chunkArray[2]=i+1;
		}
		if(array[i]=="z"){
			chunkArray[3]=i;
		}
	}
	x=Array.slice(array,chunkArray[0],chunkArray[1]);
	y=Array.slice(array,chunkArray[2],chunkArray[3]);
	setSlice(array[array.length-1]);
	makeSelection("freehand",x,y);
	run("Enlarge...", "enlarge=&enlarge");//This enlarges by scaled units which are assumed to be microns
	run("Crop");
	run("Select None");
}
/*
 * Convert array of coordinates into selections
 */
function convertArrayToRois(coordinateArray,roi){
	xList=newArray();
	yList=newArray();
	zList=newArray(); 
	for(i=0; i<coordinateArray.length; i++){
		if(coordinateArray[i]=="x"){
			xList=Array.concat(xList,i+1);
		}
		if(coordinateArray[i]=="y"){
			xList=Array.concat(xList,i);
			yList=Array.concat(yList,i+1);
		}
		if(coordinateArray[i]=="z"){
			yList=Array.concat(yList,i);
			zList=Array.concat(zList,i+1);
		}
	}
	for(i=0; i<zList.length; i++){
		index=i*2;
		xCoords=Array.slice(coordinateArray,xList[index],xList[index+1]);
		yCoords=Array.slice(coordinateArray,yList[index],yList[index+1]);
		makeSelection("Freehand", xCoords,yCoords);
		zCoords=zList[i];
		setSlice(coordinateArray[zCoords]);
		if(roi==1){
			roiManager("Add");
			run("Select None");
		}	
	}
}
/*
 * Convert ROIs into list of coordinates.
 */
function convertRoisToCoords(rois){
	selectionArray=newArray();
	if(rois==1){
		for(i=0; i<roiManager("Count"); i++){
			roiManager("Select", i);
			getSelectionCoordinates(x,y);
			z=getSliceNumber();
			selectionArray=Array.concat(selectionArray,"x",x,"y",y,"z",z);
			run("Select None");
			wait(50);
		}
		roiManager("Reset");
	}else{
		if(selectionType>0){
			getSelectionCoordinates(x,y);
			z=getSliceNumber();
			selectionArray=Array.concat(selectionArray,"x",x,"y",y,"z",z);
			run("Select None");
			roiManager("Reset");
		}
	}
	return selectionArray;
}
/*
 * Convert selections into list of coordinates.
 * Requires a binary image or a non binary image with a selection.
 */
function convertSelectionsToRois(slice,returnArray){
	arrayCoords=newArray();
	if(slice==1 || nSlices==1){
		getThreshold(min,max);
		threshold=min+max;
		if(is("binary")==true ||  threshold>-2){
			if(threshold==-2){
				setThreshold(255,255);
			}
			run("Create Selection");
			roiManager("Add");
			if(selectionType==9){
				roiManager("Split");
				roiManager("Select", 0);
				roiManager("Delete");
			}
		}
		run("Select None");
	}else{
		for(i=0; i<nSlices; i++){
			setSlice(i+1);
			getThreshold(min,max);
			threshold=min+max;
			if(is("binary")==true ||  threshold>-2){
				if(threshold==-2){
					setThreshold(255,255);
				}
				run("Create Selection");
				getStatistics(area,mean,min,max,std,histo);
				if(max>0){
					index=roiManager("Count");
					roiManager("Add");
					if(selectionType==9){
						roiManager("Split");
						roiManager("Select", index);
						roiManager("Delete");
						run("Select None");
					}
				}
				run("Select None");
			}
			resetThreshold();
		}
		resetThreshold();
	}
	run("Select None");
	resetThreshold();
	if(returnArray==1){
		for(i=0; i<roiManager("Count"); i++){
			roiManager("Select", i);
			slice=getSliceNumber();
			getSelectionBounds(x,y,width,height);
			arrayCoords=Array.concat(arrayCoords,x,y,slice);
			run("Select None");
		}
		return arrayCoords;
	}	
}
/*
 * Crop in XYZ.
 */
function cropXy(array,enlarge,allImages){
	if(allImages==0){
		run("Select None");
		convertArrayToRois(array,1);
		roiManager("Select", 0);
		run("Enlarge...", "enlarge=&enlarge");
		run("Crop");
		selection=convertRoisToCoords(0);
		run("Select None");
	}else{
		for(i=0; i<nImages; i++){
			run("Select None");
			convertArrayToRois(array,1);
			selectImage(i+1);
			roiManager("Select", 0);
			run("Enlarge...", "enlarge=&enlarge");
			run("Crop");
			selection=convertRoisToCoords(0);
			run("Select None");
		}
		for(i=0; i<nImages; i++){
			selectImage(i+1);
			run("Select None");
		}
	}
	run("Select None");

	return selection;
}
/* 
 * Return roi with greatest area.
 */
function returnLargestRoi(array,upper){
	convertArrayToRois(array,1);
	areaArray=newArray(roiManager("Count"));
	for(i=0; i<roiManager("Count"); i++){
		roiManager("Select", i);
		getStatistics(area);
		areaArray[i]=area;
		run("Select None");
	}
	Array.getStatistics(areaArray,min,max,mean,std);
	deleteArray=newArray();
	for(i=0; i<roiManager("Count"); i++){
		roiManager("Select", i);
		getStatistics(area);
		roiManager("Deselect");
		if(area!=max){
			deleteArray=Array.concat(deleteArray, i);
		}
		run("Select None");
	}
	if(deleteArray.length>0){
		roiManager("Select", deleteArray);
		roiManager("Delete");
		run("Select None");
	}
	roiManager("Select", 0);	
	getStatistics(area,mean,min,max,std,histo);	
	run("Convex Hull");
	roiManager("Update");
	roiManager("Select", 0);
	roiManager("Remove Slice Info");
	largestCoords=convertRoisToCoords(1);
	run("Select None");
	if(area<upper){
		largestCoords=newArray();
	}
	return largestCoords;
}
/*
 * Filter roiManager objects based on shape. 
 */
function shapeFilterRois(array,areaferet,feretcirc){
	convertArrayToRois(array,1);
	if(roiManager("Count")>0){
		deleteArray=newArray();
		for(i=0; i<roiManager("Count"); i++){
			roiManager("Select", i);
			getStatistics(area);
			List.setMeasurements();
			feret=List.getValue("Feret");
			circ=List.getValue("Circ.");
			if(areaferet==1){
				if(area>2000 || feret>60){
					deleteArray=Array.concat(deleteArray,i);
				}
			}
			if(feretcirc==1){
				if(feret>47 && circ<0.85){
					if(deleteArray.length>0){
						if(deleteArray[deleteArray.length-1]!=i){
							deleteArray=Array.concat(deleteArray, i);
						}
					} else{
					deleteArray=Array.concat(deleteArray, i);
					}
				}
			}
			run("Select None");
		}
		if(deleteArray.length>0){
			roiManager("Select", deleteArray);
			roiManager("Delete");
			run("Select None");
		}
	}
	if(roiManager("Count")>0){
		array=convertRoisToCoords(1);
	}
	return array;
}
/*
 * A general method to extract features from images.
 * Gaussian is the size of the blurring radius in microns (3 micron works well to find Tet)
 * Laplacian is the size of the Laplacian radius in microns (1
 */
function roughFeatureExtraction(gaussian, laplacian){
	getVoxelSize(width,height,depth,unit);
	run("Z Project...", "projection=[Max Intensity]");
	temp=getImageID;
	getStatistics(area,mean,min,max,std);
	run("Subtract...", "value=&mean");
	run("Gaussian Blur...", "sigma=&gaussian scaled");
	factor=(0.12578/width)*laplacian;
	run("FeatureJ Laplacian", "compute smoothing=&factor");
	temp2=getImageID;
	run("Invert");
	run("Select All");
	shrink=-1*gaussian/2;
	run("Enlarge...", "enlarge=&shrink");
	getStatistics(area,mean,min,max,std,histogram);
	setMinAndMax(min,max);
	run("Select None");
	run("8-bit");
	run("Select All");
	run("Enlarge...", "enlarge=&shrink");
	run("Gaussian Blur...", "sigma=&gaussian scaled");
	setAutoThreshold("Triangle dark");
	getThreshold(threshold,max);
	getStatistics(area,mean,min,max,std,histogram);
	setThreshold(threshold, 255);
	if(max>=threshold){
		convertSelectionsToRois(1,0);
		features=convertRoisToCoords(1);
	}
	selectImage(temp);
	close();
	selectImage(temp2);
	close();

	return features;
}
/*
 * Determine if image has a selection. If image does not have a selection
 * it returns 1. If image has a selection it returns 0.
 */
function queryForImageSelection(){
	getStatistics(area,mean,min,max,std,histo);
	getVoxelSize(width,height,depth,unit);
	wide=getWidth()*width;
	high=getHeight()*height;
	area2=wide*high;
	if(round(area)==round(area2)){
		return 1;
	} else{
		return 0;
	}
}
/*
 * Find all local maxima within a 3D stack
 */
function findLocalMaxima(selection,constant){
	dup="Duplicate";
	run("Duplicate...", "title=[&dup] duplicate");
	dup=getImageID();
	dup2="Duplicate2";
	run("Duplicate...", "title=[&dup2] duplicate");
	dup2=getImageID();
	getVoxelSize(height,width,depth,unit);
	blurRadius=0.125/width;
	if(blurRadius<1){
		blurRadius=1;
	}
	//run("Gaussian Blur...", "sigma=&blurRadius stack");
	maxXyRadius=0.5/width;
	if(maxXyRadius<1){
		maxXyRadius=1;
	}
	maxZRadius=1.2/depth;
	if(maxZRadius<1){
		maxZRadius=1;
	}
	if(selection.length>0){
		convertArrayToRois(selection,0);
		run("Enlarge...", "enlarge=-1 pixel");
		run("Clear Outside","stack");
		run("Select None");
	}
	run("3D Fast Filters","filter=MaximumLocal radius_x_pix=&maxXyRadius radius_y_pix=&maxXyRadius radius_z_pix=&maxZRadius Nb_cpus=4");
	temp=getImageID();
	noSlices=nSlices;
	makeBinaryImage();
	convertSelectionsToRois(0,0);
	/*
	 * first round
	 */
	selectImage(dup);
	intensity=getRoiIntensity();
	Array.getStatistics(intensity,min,max,globalmean,globalstd);
	globalThreshold=globalmean+globalstd;
	slices=getRoiSlice();
	rollingAvgArray=rollingSliceAverage(5,noSlices);
	thresholdArray=getRoiAvgWithinRange(intensity,slices,rollingAvgArray,noSlices,constant);
		for(i=0; i<thresholdArray.length; i++){
		setResult("threshold", i, thresholdArray[i]);
	}
	updateResults;
	deleteRoisBelowThreshold(intensity,thresholdArray,globalThreshold);
	consolidatMultiVoxelObjects();
	selectImage(dup);
	coords=convertRoisToCoords(1);	
	selectImage(temp);
	close();
	selectImage(dup);
	close();
	selectImage(dup2);
	close();
	roiManager("Reset");

	return coords;
}
/*
 * Convert Roi to Binary Image.
 * Requires that an appropriately sized image is selected.
 * It returns a binary image title and image
 */
function convertRoiToBinary(){
	run("Select None");
	dup="Binary_"+getTitle();;
	run("Duplicate...", "title=[&dup] duplicate");
	dupID=getImageID();
	run("Select All");
	run("Clear","stack");
	run("Select None");
	run("8-bit");
	for(i=0; i<roiManager("Count"); i++){
		roiManager("Select", i);
		setColor(255);
		fill();
		run("Select None");
	}
	makeBinaryImage();
	return dupID;
}
/*
 * Consolidate multi-voxel objects to centroids.
 */
function consolidatMultiVoxelObjects(){
	binaryID=convertRoiToBinary();
	run("3D OC Options", "  dots_size=1 font_size=10 redirect_to=none");
	run("3D Objects Counter", "threshold=128 slice=10 min.=1 max.=20000000 centroids");
	consolidatedObjects=getImageID();
	makeBinaryImage();
	roiManager("Reset");
	convertSelectionsToRois(0,0);
	selectImage(consolidatedObjects);
	close();
	selectImage(binaryID);
	close();
}
/*
 * Make any image a binary image.
 */
function makeBinaryImage(){
	if(bitDepth==16){
		setMinAndMax(0,1);
		run("8-bit");
		setMinAndMax(0,0);
		run("Apply LUT", "stack");
	}
	if(bitDepth==8){
		setMinAndMax(0,0);
		run("Apply LUT", "stack");
	}
}
/*
 * This loop will get the maximum intensity for every object
 */
function getRoiIntensity(){
	objects=roiManager("Count");
	objectIntensityArray=newArray(objects);
	for(i=0; i<objects; i++){
		roiManager("Select", i);
		getStatistics(area,mean,min,max,std,histo);
		objectIntensityArray[i]=mean;
	}
	return objectIntensityArray;
}
/*
 * This loop will get the slice for every object
 */
function getRoiSlice(){
	objects=roiManager("Count");
	objectSliceArray=newArray(objects);
	for(i=0; i<objects; i++){
		roiManager("Select", i);
		slice=getSliceNumber();
		objectSliceArray[i]=slice;
	}
	return objectSliceArray;
}
/*
 * This loop will create an array that contains the start and stop slice
 * for going through a stack and averaging things over a 5 slice interval
 * centered around the current slice.
 */ 
function rollingSliceAverage(range,slices){
	rollingAvgArray=newArray(slices*2);
	tempArray=newArray(range*2+1);
	counter=-1*range;
	for(ii=0; ii<range*2+1; ii++){
		tempArray[ii]=counter;
		counter++;
	}
	counter=0;
	for(i=1; i<slices+1; i++){
		tempArray2=newArray();
		for(ii=0; ii<tempArray.length; ii++){
			if(tempArray[ii]+i>0 && tempArray[ii]+i<slices+1){
				value=tempArray[ii]+i;
				tempArray2=Array.concat(tempArray2,value);
			}
		}
		Array.getStatistics(tempArray2,min,max,mean,std);
		rollingAvgArray[counter]=min;
		counter++;
		rollingAvgArray[counter]=max;
		counter++;
	}
	return rollingAvgArray;
}
/*
 * This loop will go through and get the average and standard deviation
 * for all objects within the slice range specified usng the loop above.
 * It will return an object threshold for each slice that is tailored to the
 * objects within that slice and the objects within the rolling slice average.
 * array1 is the array of object intensities.
 * array2 is the array of object slices.
 * array3 is the array containing the start and stop for rolling average
 */
function getRoiAvgWithinRange(array1,array2,array3,slices,constant){
	objects=array1.length;
	thresholdArray=newArray(slices);
	counter=0;
	for(i=0; i<slices; i++){
		index=i*2;
		tempArray=newArray(objects);
		counter1=0;
		counter2=0;
		count1Reset=0;
		for(ii=0; ii<objects; ii++){
			if(array2[ii]>=array3[index] && array2[ii]<=array3[index+1]){
				if(count1Reset==0){
					counter1=ii;
					count1Reset=1;
				}
				tempArray[ii]=array1[ii];
			} else if(array2[ii]>array3[index+1]){	
				counter2=ii;
				ii=objects;
			}
		}
		if(array3[index+1]==slices){
			counter2=array2.length-1;
		}
		tempArray2=Array.slice(tempArray,counter1,counter2);
		Array.getStatistics(tempArray2,min,max,mean,std);
		cv=std/mean;//This needs to be refined for EM gain images
		thresholdArray[counter]=((constant-(cv*constant))*std)+mean;
		if(cv>1){
			thresholdArray[counter]=std+mean;
		}
		counter++;
	}
	return thresholdArray;
}
/*
 * This loop will go through and delete all ROI that are below
 * their slice range threshold.
 * array1 is a list of object intensities
 * array2 is a list of slice thresholds
 */
function deleteRoisBelowThreshold(array1,array2,singleThreshold){
	objects=roiManager("Count");
	coordArray=newArray();
	
	for(i=0; i<objects; i++){
		roiManager("Select",i);
		getSelectionCoordinates(x,y);
		slice=getSliceNumber();
		threshold=array2[slice-1];
		threshold2=singleThreshold;
		if(array1[i]>threshold && array1[i]>threshold2){
			coordArray=Array.concat(coordArray,slice);
			coordArray=Array.concat(coordArray,x[0]);
			coordArray=Array.concat(coordArray,y[0]);
		}
	}
	roiManager("Reset");
	for(i=0; i<coordArray.length/3; i++){
		index=i*3;
		setSlice(coordArray[index]);
		makeRectangle(coordArray[index+1]-1,coordArray[index+2]-1,1,1);
		roiManager("Add");
	}
}
/*
 * Remove x,y,z from an array of coordinates.
 * This is only useful if the coordinates are points because the spacing between
 * objects is assumed to be 1
 */
function cleanUpCoordArray(array){
	cleanArray=newArray(); 
	for(i=0; i<array.length; i++){
		if(array[i]=="x"){
			cleanArray=Array.concat(cleanArray,parseFloat(array[i+1]));
		}
		if(array[i]=="y"){
			cleanArray=Array.concat(cleanArray,parseFloat(array[i+1]));
		}
		if(array[i]=="z"){
			cleanArray=Array.concat(cleanArray,parseFloat(array[i+1]));
		}
	}
	return cleanArray;
}
/*
 * Remove BBs from the middle. A convex hull image from the average BBs will be created and a 3D EDM of this image will be created.
 */
function removeMidBbs(distance){
	getVoxelSize(width,height,depth,unit);
	dup="Duplicate";
	run("Duplicate...", "title=[&dup] duplicate");
	dup=getImageID();
	run("Select All");
	run("Clear","stack");
	run("8-bit");
	convertRoiToBinary();
	binary=getImageID();
	setSlice(1);
	run("Select All");
	getStatistics(area,mean,min,max,std,histo);
	if(max>0){
		run("Clear","slice");
	}
	run("Select None");
	waitForUser("1-1");
	run("3D Convex Hull");
	convex=getImageID();
	
	run("3D Distance Map", "threshold=1");
	distanceMap=getImageID();
	temp=newArray();
	for(i=0; i<roiManager("Count"); i++){
		roiManager("Select", i);
		getStatistics(area,mean,min,max,std,histo);
		if(mean>distance){
			temp=Array.concat(temp,i);
		}
		run("Select None");
	}
	if(temp.length>0){
		roiManager("Select", temp);
		roiManager("Delete");	
	}
	selectWindow("Log");
	run("Close");
	selectImage(dup);
	close();
	selectImage(binary);
	close();
	selectImage(convex);
	close();
	selectImage(distanceMap);
	close();
}
/*
 * Filter maxima based on skew and kurtosis. Also returns the 3D coordinate of the OA.
 * The BB image must be selected heading into the macro.
 */
function skewKurtFilterHomog(enlargeBox){
	getVoxelSize(width,height,depth,unit);
	title=getTitle();
	dup="Duplicate";
	run("Duplicate...", "title=[&dup] duplicate");
	dup=getImageID();
	for(i=0; i<nSlices; i++){
		threshold=getResult("threshold",i);
		setSlice(i+1);
		run("Subtract...", "value=&threshold slice");
	}
	makeBinaryImage();
	minSize=0.5/(width*height*depth);
	run("3D OC Options", "volume integrated_density centre_of_mass dots_size=1 font_size=10 redirect_to=[&title]");
	run("3D Objects Counter", "threshold=128 slice=26 min.=&minSize max.=4530864 objects statistics");
	oaObjects=getImageID();
	selectImage(dup);
	close();
	intDenArray=newArray(nResults);
	for(i=0; i<nResults; i++){
		intDenArray[i]=getResult("IntDen",i);
	}
	Array.getStatistics(intDenArray,min,max,mean,std);
	index=100000000000;
	for(i=0; i<nResults; i++){
		value=intDenArray[i];
		if(value==max){
			index=i;
		}
	}
	oaX=getResult("XM",index);
	oaY=getResult("YM",index);
	oaZ=getResult("ZM",index);
	oaArray=newArray(oaX,oaY,oaZ);
	if(isOpen("Results")==1){
		selectWindow("Results");
		run("Close");
	}
	deleteArray=newArray();
	for(i=0; i<roiManager("Count"); i++){
		selectImage(title);
		roiManager("Select", i);
		slice=getSliceNumber();
		factor=round((0.1278/width)*enlargeBox);
		run("Enlarge...", "enlarge=&factor pixel");
		run("To Bounding Box");
		List.setMeasurements();
		skew=List.getValue("Skew");
		kurt=List.getValue("Kurt");
		if(skew+kurt<1){
			selectImage(oaObjects);
			roiManager("Select", i);
			getStatistics(area,mean,min,max,mean,std);
			if(min>0){
				deleteArray=Array.concat(deleteArray, i);
			}
		}
		run("Select None");
	}
	selectImage(oaObjects);
	close();
	if(deleteArray.length>0){
		roiManager("Select", deleteArray);
		roiManager("Delete");
		roiManager("Deselect");
		run("Select None");
	}
	coords=convertRoisToCoords(1);
	convertArrayToRois(coords,1);
	returnCoords=cleanUpCoordArray(coords);
	returnArray=Array.concat(returnCoords,oaArray);
	return returnArray;
}
/*
 * This plots BB average location
 */
function plotAverageBBLocation(objectArray,neighborArray,cluster,plot){ 
	getVoxelSize(width,height,depth,unit);
	if(plot==1){
		dup="Duplicate";
		run("Duplicate...", "title=[&dup] duplicate");
		dup=getImageID();
		run("Select All");
		run("Clear","stack");
		run("8-bit");
	}
	averageArray=newArray(objectArray.length);
	for(i=0; i<objectArray.length/3; i++){
		index=i*3;
		x1=objectArray[index];
		y1=objectArray[index+1];
		z1=objectArray[index+2];
		index2=i*cluster*3;
		x2=0;
		y2=0;
		z2=0;
		for(ii=0; ii<cluster; ii++){
			index3=index2+(ii*3);
			x2+=neighborArray[index3];
			y2+=neighborArray[index3+1];
			z2+=neighborArray[index3+2];
		}
		x=x2/cluster;
		y=y2/cluster;
		z=z2/cluster;
		averageArray[index]=x;
		averageArray[index+1]=y;
		averageArray[index+2]=z;
		if(plot==1){
			setSlice(z+1);
			makeRectangle(x,y,1,1);
			setColor(255);
			fill();
			run("Select None");
		}
	}
	if(plot==1){
		selectImage(dup);
		close();
	}
	return averageArray;
}
/*
 * Find most distant partner for each BB
 */
function findPoles(cluster,array,group){
	getVoxelSize(width,height,depth,unit);
	clusterCoordinatesArray=newArray();//This will be the finall array that has all the coordinates for the clustered objects
	maxDistanceArray=newArray();
	for(i=0; i<array.length/3; i++){//This marches through all the objects of interest
		index=i*3;//0,3,6,9,12...This is the index for pulling the coordinates for the object from the coordinatesArray
		x1=array[index]*width;//0,3,6,9,12...x coordinate for the object
		y1=array[index+1]*height;//1,4,7,10,13...y coordinate for the object
		z1=array[index+2]*depth;//2,5,8,11,14...z coordinate for the object
		distanceArray=newArray(array.length/3);//This array contains all the distances relative to the object of interest
		for(ii=0; ii<array.length/3; ii++){//This marches through all objects a second time
			index2=ii*3;//0,3,6,9,12...This index is the index for pulling the objects the coordinates for the object from the coordinatesArray
			x2=array[index2]*width;//0,3,6,9,12...x coordinate for the 2nd object
			y2=array[index2+1]*height;//1,4,7,10,13...y coordinate for the 2nd object
			z2=array[index2+2]*depth;//2,5,8,11,14...z coordinate for the 2nd object
			distance=sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));//This is the distance between the first and second object
			distanceArray[ii]=distance;//This assigns the distance between the two objects to the distanceArray
		}
		ranked=Array.rankPositions(distanceArray);//This ranks the distanceArray so that a new array is created that contains the relative positions of the 2nd objects relative to the first object
		rankedDistanceArray=Array.rankPositions(ranked);
		tempClusterCoordinatesArray=newArray(cluster*3);//This creates a new array that has enough space to hold the xyz coordinates for nth closest objects
		tempMaxDistanceArray=newArray(cluster);
		counter=0;
		for(ii=0; ii<array.length/3; ii++){//This loop marches through the rankedDistanceArray looking for positions that are less than cluster size
			value=rankedDistanceArray[ii];
			if(value==rankedDistanceArray.length-1){//If the rankedDistance is less than the cluster size then extract the spatial coordinates
				index2=counter*3;
				tempClusterCoordinatesArray[index2]=array[(ii*3)];
				tempClusterCoordinatesArray[index2+1]=array[(ii*3)+1];
				tempClusterCoordinatesArray[index2+2]=array[(ii*3)+2];
				tempMaxDistanceArray[counter]=distanceArray[ii];
				counter++;
			}
		}
		clusterCoordinatesArray=Array.concat(clusterCoordinatesArray,tempClusterCoordinatesArray);
		maxDistanceArray=Array.concat(maxDistanceArray,tempMaxDistanceArray);
	}
	ranked=Array.rankPositions(maxDistanceArray);
	rankedMaxDistanceArray=Array.rankPositions(ranked);
	counter=0;
	avgSize=group;
	poleArray=newArray((avgSize)*6);
	xArray=newArray((avgSize)*2);
	for(i=0; i<rankedMaxDistanceArray.length; i++){
		if(rankedMaxDistanceArray[i]>rankedMaxDistanceArray.length-(avgSize+1)){				
			index=counter*6;
			poleArray[index]=array[i*3];
			poleArray[index+1]=array[(i*3)+1];
			poleArray[index+2]=array[(i*3)+2];
			poleArray[index+3]=clusterCoordinatesArray[i*3];
			poleArray[index+4]=clusterCoordinatesArray[(i*3)+1];
			poleArray[index+5]=clusterCoordinatesArray[(i*3)+2];
			counter++;
		}
	}
	Array.getStatistics(xArray,min,max,mean,std);
	pole1X=0;
	pole1Y=0;
	pole1Z=0;
	pole2X=0;
	pole2Y=0;
	pole2Z=0;
	counter1=0;
	counter2=0;
	counter=0;
	poleStartX=poleArray[0]*width;
	poleStartY=poleArray[1]*height;
	poleStartZ=poleArray[2]*depth;
	for(i=0; i<poleArray.length/3; i++){
		index=i*3;
		distance=sqrt(pow(poleArray[index]*width-poleStartX,2)+pow(poleArray[index+1]*height-poleStartY,2)+pow(poleArray[index+2]*depth-poleStartZ,2));
		if(distance<15){
			pole1X+=poleArray[index];
			pole1Y+=poleArray[index+1];
			pole1Z+=poleArray[index+2];
			counter1++;
		} else{
			pole2X+=poleArray[index];
			pole2Y+=poleArray[index+1];
			pole2Z+=poleArray[index+2];
			counter2++;
		}
		counter++;
	}
	pole1X=pole1X/counter1;
	pole1Y=pole1Y/counter1;
	pole1Z=pole1Z/counter1;
	pole2X=pole2X/counter2;
	pole2Y=pole2Y/counter2;
	pole2Z=pole2Z/counter2;
	run("Specify...", "width=30 height=30 x=&pole1X y=&pole1Y slice=&pole1Z centered");
	getStatistics(area,mean1,min,max,std);
	run("Select None");
	run("Specify...", "width=30 height=30 x=&pole2X y=&pole2Y slice=&pole2Z centered");
	getStatistics(area,mean2,min,max,std);
	run("Select None");
	poleArray=newArray(6);
	if(mean1>mean2){
		poleArray[0]=pole1X;
		poleArray[1]=pole1Y;
		poleArray[2]=pole1Z;
		poleArray[3]=pole2X;
		poleArray[4]=pole2Y;
		poleArray[5]=pole2Z;
	} else{
		poleArray[0]=pole2X;
		poleArray[1]=pole2Y;
		poleArray[2]=pole2Z;
		poleArray[3]=pole1X;
		poleArray[4]=pole1Y;
		poleArray[5]=pole1Z;
	}
	return poleArray;
}
function crossProduct(array){
	vecA1=array[0]-(array[3]);
	vecA2=array[1]-(array[4]);
	vecA3=array[2]-(array[5]);
	vecB1=array[0]-(array[6]);
	vecB2=array[1]-(array[7]);
	vecB3=array[2]-(array[8]);
	vecC1=vecA2*vecB3-vecA3*vecB2;
	vecC2=vecA3*vecB1-vecA1*vecB3;
	vecC3=vecA1*vecB2-vecA2*vecB1;
	cp=newArray(vecC1,vecC2,vecC3);
	return cp;
}
function createPlaneEquation(array){
	xConstant=array[0];
	yConstant=array[1];
	zConstant=array[2];
	offset=(array[0]*(-1*array[3]))+(array[1]*(-1*array[4]))+(array[2]*(-1*array[5]));
	planeEquation=newArray(xConstant,yConstant,zConstant,offset);
	return planeEquation;
}
function distPointPlane(array){
	topTerm=abs((array[3]*array[0])+(array[4]*array[1])+(array[5]*array[2])+(array[6]));
	bottomTerm=sqrt(pow(array[3],2)+pow(array[4],2)+pow(array[5],2));
	distance=topTerm/bottomTerm;
	return distance;
}
function planePenalty(objects,neighbors,poles){
	getVoxelSize(width,height,depth,unit);
	clusterSize=neighbors.length/objects.length;
	arrayReturn=newArray();
	counter=0;
	for(i=0; i<objects.length/3; i++){
		index=i*3;
		point=newArray(objects[index],objects[index+1],objects[index+2]);
		array=Array.concat(poles,point);
		cp=crossProduct(array);
		array=Array.concat(cp,point);
		equation=createPlaneEquation(array);
		tempArrayReturn=newArray(clusterSize);
		index2=i*clusterSize*3;
		for(ii=0; ii<clusterSize; ii++){
			index3=index2+(ii*3);
			point2=newArray(neighbors[index3],neighbors[index3+1],neighbors[index3+2]);
			array=Array.concat(point2,equation);
			pointPlaneDistance=distPointPlane(array);
			//arrayReturn[counter]=pointPlaneDistance;
			tempArrayReturn[ii]=pointPlaneDistance;
			counter++;
		}
		rank=Array.rankPositions(tempArrayReturn);
		rank2=Array.rankPositions(rank);
		arrayReturn=Array.concat(arrayReturn,rank2);
	}
	return arrayReturn;
}
function partnerDistancePenalty(objects,neighbors){
	getVoxelSize(width,height,depth,unit);
	clusterSize=neighbors.length/objects.length;
	arrayReturn=newArray((objects.length/3)*(clusterSize));
	counter=0;
	for(i=0; i<objects.length/3; i++){
		index=i*3;
		x1=objects[index]*width;//0,3,6,9,12...x coordinate for the object
		y1=objects[index+1]*height;//1,4,7,10,13...y coordinate for the object
		z1=objects[index+2]*depth;//2,5,8,11,14...z coordinate for the object
		index2=i*clusterSize*3;
		for(ii=0; ii<clusterSize; ii++){
			index3=index2+(ii*3);
			x2=neighbors[index3]*width;
			y2=neighbors[index3+1]*height;
			z2=neighbors[index3+2]*depth;
			distance=sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));
			arrayReturn[counter]=distance;
			if(distance==0){
				arrayReturn[counter]=10000;
			}
			counter++;
		}
	}
	return arrayReturn;
}
function polePenalty(objects,neighbors,poles,ant){
	getVoxelSize(width,height,depth,unit);
	clusterSize=neighbors.length/objects.length;
	arrayReturn=newArray((objects.length/3)*(clusterSize));
	if(ant==1){
		antx=(poles[0]*width)+(-1*((poles[3]*width)-(poles[0]*width)));
		anty=(poles[1]*width)+(-1*((poles[4]*width)-(poles[1]*width)));
		antz=(poles[2]*depth);
	} else{
		antx=(poles[3]*width)+(-1*((poles[0]*width)-(poles[3]*width)));
		anty=(poles[4]*width)+(-1*((poles[1]*width)-(poles[4]*width)));
		antz=(poles[5]*depth);
	}
	counter=0;
	for(i=0; i<objects.length/3; i++){
		index=i*3;
		x1=objects[index]*width;
		y1=objects[index+1]*height;
		z1=objects[index+2]*depth;
		distance1=sqrt(pow(antx-x1,2)+pow(anty-y1,2)+pow(antz-z1,2));
		index2=i*clusterSize*3;
		for(ii=0; ii<clusterSize; ii++){
			index3=index2+(ii*3);
			x2=neighbors[index3]*width;
			y2=neighbors[index3+1]*height;
			z2=neighbors[index3+2]*depth;			
			distance2=sqrt(pow(antx-x2,2)+pow(anty-y2,2)+pow(antz-z2,2));
			if(distance1>distance2){
				arrayReturn[counter]=1;
			} else{
				arrayReturn[counter]=10000;
			}
			counter++;
		}
	}
	return arrayReturn;
}
/*
 * Sort scores to find the lowest score. 
 * Return the components of the lowest score along with the xyz coordinates for the lowest score
 */
function sortScores(cluster,score1,score2,score3,neighbors,points){
	getVoxelSize(width,height,depth,unit);
	arrayReturn=newArray((score1.length/cluster)*6);//This array will return the XYZ and individual score components for every object along with a binary value indicating whether the object 
	scoreArray=newArray(score1.length);
	for(i=0; i<score1.length; i++){
		scoreArray[i]=score1[i]*score2[i]*score3[i];
	}
	for(i=0; i<scoreArray.length/cluster; i++){//This loop goes through and connects everything based solely on its score
		index=i*cluster; 
		index2=i*6;
		index3=index*3;
		tempArray=newArray(cluster);
		for(ii=0; ii<cluster; ii++){
			tempArray[ii]=scoreArray[index+ii];
		}
		ranked=Array.rankPositions(tempArray);
		rankedScores=Array.rankPositions(ranked);
		for(ii=0; ii<cluster; ii++){
			if(rankedScores[ii]==1){
				arrayReturn[index2]=neighbors[((index+ii)*3)];
				arrayReturn[index2+1]=neighbors[((index+ii)*3)+1];
				arrayReturn[index2+2]=neighbors[((index+ii)*3)+2];
				arrayReturn[index2+3]=score1[index+ii];
				arrayReturn[index2+4]=score2[index+ii];
				arrayReturn[index2+5]=score3[index+ii];
			}
		}
	}
	deletePositions=newArray();
	for(i=0; i<arrayReturn.length/6; i++){//This loop goes back and finds anything that's received two connections
		index=i*6;
		index2=i*3;
		objectX=points[index2];
		objectY=points[index2+1];
		objectZ=points[index2+2];
		tempX=arrayReturn[index];//X coordinate of the object
		tempY=arrayReturn[index+1];//Y coordinate of the object
		tempZ=arrayReturn[index+2];//Z coordinate of the object
		tempScore=arrayReturn[index+4];//Distance between the object and its connected partner which is used to break ties
		duplicate=0;//This will stay 0 unless the search finds the same object more than once
		tempScoreSort=newArray();//An array to hold the scores of objects that are connected more than once
		tempPositionSort=newArray();//An array to hold the positions of the objects that are connected more than once. This should be the index
		for(ii=0; ii<arrayReturn.length/6; ii++){//This loop starts searching the partners for connections that were made more than once.
			index3=ii*6;
			tempX1=arrayReturn[index3];//X coordinate of the connected object
			tempY1=arrayReturn[index3+1];//Y coordinate of the connected object
			tempZ1=arrayReturn[index3+2];//Z coordinate of the connected object
			tempScore1=arrayReturn[index3+4];//Distance between the connected object and its partner
			if(ii!=i){//This prevents the search from finding the connected object itself
				if(tempX==tempX1 && tempY==tempY1 && tempZ==tempZ1){//This is only satisfied if it found the same object
					duplicate=1;//This lets the loop know that it found the same connected object
					tempScoreSort=Array.concat(tempScoreSort,tempScore1);//This adds the score to the score array
					tempPositionSort=Array.concat(tempPositionSort,index3);//This adds the position of multiple connected objects to the position array
				}
			}
		}
		if(duplicate==1){//If a multiple connected object was found than add the initial object to the array
			tempScoreSort=Array.concat(tempScoreSort,tempScore);
			tempPositionSort=Array.concat(tempPositionSort,index);
			rank=Array.rankPositions(tempScoreSort);
			rankedScores=Array.rankPositions(rank);
			for(ii=0; ii<rankedScores.length; ii++){
				if(rankedScores[ii]!=0){
					deleteIndex=tempPositionSort[ii];
					deletePositions=Array.concat(deletePositions,deleteIndex);
				}
			}
		}
		duplicate=0;
	}
	for(i=0; i<deletePositions.length; i++){
		index=deletePositions[i];
		arrayReturn[index]=-1;
		arrayReturn[index+1]=-1;
		arrayReturn[index+2]=-1;
		arrayReturn[index+3]=-1;
		arrayReturn[index+4]=-1;
		arrayReturn[index+5]=-1;
	}
	return arrayReturn;
}
/*
 * Find points to reconnect
 */
function findPointsToReconnect(connected,points){
	getVoxelSize(width,height,depth,unit);
	reconnectArray=newArray();
	for(i=0; i<connected.length/6; i++){//This loop goes through and connects everything based solely on its score
		index=i*6;
		index2=i*3;
		objectx=points[index2];
		objecty=points[index2+1];
		objectz=points[index2+2];		
		if(connected[index]==-1){
			reconnectArray=Array.concat(reconnectArray,objectx,objecty,objectz);
		}
	}
	return reconnectArray;
}
/*
 * Find possible partners to reconnect to
 */
function findPartnersToReconnect(connected,points){
	getVoxelSize(width,height,depth,unit);
	possibleArray=newArray();
	for(i=0; i<connected.length/6; i++){//This loop goes through and connects everything based solely on its score
		index=i*6;
		index2=i*3;
		objectx=points[index2];
		objecty=points[index2+1];
		objectz=points[index2+2];
		possible=0;
		for(ii=0; ii<connected.length/6; ii++){
			index3=ii*6;
			queryx=connected[index3];
			queryy=connected[index3+1];
			queryz=connected[index3+2];
			if(objectx==queryx && objecty==queryy && objectz==queryz){
				possible=1;
			}
		}
		if(possible==0){
			possibleArray=Array.concat(possibleArray,objectx,objecty,objectz);
		}
	}
	return possibleArray;
}
/*
 * Iterative basal body finding algorithm with using planes
 */
function iterativeBasalBodyFinding(points,clusterSize,poles){
	row=0;
	scoring=1;
	intensity=1;
	totalObjects=points.length/3;
	finalConnections=newArray();
	for(i=0; i<100; i++){
		startRow=row;
		if(i==0){
			allConnections=makeConnectionsShell(points,15,poles,scoring);
		} else{
			allConnections=makeConnectionsShellLoop(antReconnectPoints,postReconnectPoints,poles,scoring);
		}
		anteriorConnections=Array.slice(allConnections,0,allConnections.length/2);
		posteriorConnections=Array.slice(allConnections,allConnections.length/2,allConnections.length);
		robustConnectionsAnt=keepRobustConnections(anteriorConnections,posteriorConnections);
		robustConnectionsPost=keepRobustConnections(posteriorConnections,anteriorConnections);
		antReconnectPoints=createNewPointList(robustConnectionsAnt,anteriorConnections);
		postReconnectPoints=createNewPointList(robustConnectionsPost,posteriorConnections);
		for(ii=0; ii<robustConnectionsAnt.length; ii++){
			if(robustConnectionsAnt[ii]==1){
				index=ii*9;				
				setResult("X",row,anteriorConnections[index]);
				setResult("Y",row,anteriorConnections[index+1]);
				setResult("Z",row,anteriorConnections[index+2]);
				setResult("X1",row,anteriorConnections[index+3]);
				setResult("Y1",row,anteriorConnections[index+4]);
				setResult("Z1",row,anteriorConnections[index+5]);				
				finalConnections=Array.concat(finalConnections,anteriorConnections[index],anteriorConnections[index+1],anteriorConnections[index+2],anteriorConnections[index+3],anteriorConnections[index+4],anteriorConnections[index+5]);
				row++;
			}
		}
		stopRow=row;
		if(startRow==stopRow){
			scoring++;
			i=100;
			antReconnectPoints1=createNewPointList(robustConnectionsAnt,anteriorConnections);
			postReconnectPoints1=createNewPointList(robustConnectionsPost,posteriorConnections);
		}
	}
	for(i=0; i<100; i++){
		startRow=row;
		if(i==0){
			finishConnectionsArray=finishConnections(antReconnectPoints1,postReconnectPoints1,poles);
		}else{
			finishConnectionsArray=finishConnections(antReconnectPoints2,postReconnectPoints2,poles);	
		}
		cleanUpConnection=cleanUpConnections(finishConnectionsArray);
		for(ii=0; ii<cleanUpConnection.length/6; ii++){
			index=ii*6;
			setResult("X",row,cleanUpConnection[index]);
			setResult("Y",row,cleanUpConnection[index+1]);
			setResult("Z",row,cleanUpConnection[index+2]);
			setResult("X1",row,cleanUpConnection[index+3]);
			setResult("Y1",row,cleanUpConnection[index+4]);
			setResult("Z1",row,cleanUpConnection[index+5]);
			finalConnections=Array.concat(finalConnections,cleanUpConnection[index],cleanUpConnection[index+1],cleanUpConnection[index+2],cleanUpConnection[index+3],cleanUpConnection[index+4],cleanUpConnection[index+5]);
			row++;
		}
		antReconnectPoints2=createNewPointList2(points,finalConnections,1);
		postReconnectPoints2=createNewPointList2(points,finalConnections,0);
		stopRow=row;
		if(startRow==stopRow){
			i=100;
		}
	}
	return finalConnections;
}

/*
 * create new point list 2
 */
function createNewPointList2(array,array2,antpost){
	reconnect=newArray();
	if(antpost==1){
		for(i=0; i<array.length/3; i++){
			index=i*3;
			x=array[index];
			y=array[index+1];
			z=array[index+2];
			found=0;
			for(ii=0; ii<array2.length/6; ii++){
				index2=ii*6;
				x2=array2[index2];
				y2=array2[index2+1];
				z2=array2[index2+2];
				if(x==x2&&y==y2&&z==z2){
					found=1;
				}
			}
			if(found==0){
				reconnect=Array.concat(reconnect,x,y,z);
			}
		}
	}
	if(antpost==0){
		for(i=0; i<array.length/3; i++){
			index=i*3;
			x=array[index];
			y=array[index+1];
			z=array[index+2];
			found=0;
			for(ii=0; ii<array2.length/6; ii++){
				index2=ii*6;
				x2=array2[index2+3];
				y2=array2[index2+4];
				z2=array2[index2+5];
				if(x==x2&&y==y2&&z==z2){
					found=1;
				}
			}
			if(found==0){
				reconnect=Array.concat(reconnect,x,y,z);
			}
		}
	}
	return reconnect;
}
/*
 * Finish connections. Array 1 is things to connect. Array 2 is things to connect to.
 */
function finishConnections(array1,array2,poles){
	getVoxelSize(width,height,depth,unit);
	antx=poles[0];
	anty=poles[1];
	antz=poles[2];
	connections=newArray();
	returnArray=newArray();
	for(i=0; i<array1.length/3; i++){
		index=i*3;
		x=array1[index];
		y=array1[index+1];
		z=array1[index+2];
		poleDist=sqrt(pow(x*width-antx*width,2)+pow(y*height-anty*height,2)+pow(z*depth-antz*depth,2));
		tempArray=newArray(array2.length/3);
		for(ii=0; ii<array2.length/3; ii++){
			index2=ii*3;
			x2=array2[index2];
			y2=array2[index2+1];
			z2=array2[index2+2];
			distance=sqrt(pow(x*width-x2*width,2)+pow(y*height-y2*height,2)+pow(z*depth-z2*depth,2));
			poleDist2=sqrt(pow(x2*width-antx*width,2)+pow(y2*height-anty*height,2)+pow(z2*depth-antz*depth,2));
			tempArray[ii]=distance;
			if(distance==0){
				tempArray[ii]=10000000000000;
			}
			if(poleDist2>poleDist){
				tempArray[ii]=10000000000000;
			}
		}
		Array.getStatistics(tempArray,min,max,mean,std);
		for(ii=0; ii<tempArray.length; ii++){
			value=tempArray[ii];
			if(value==min && value<5){
				connections=Array.concat(connections,x,y,z,array2[(ii*3)],array2[(ii*3)+1],array2[(ii*3)+2],min);
			}
		}
	}	
	return connections;
}
/*
 * Clean up connections
 */
function cleanUpConnections(array){
	getVoxelSize(width,height,depth,unit);
	cleanedUpConnections=newArray();
	usedArray=newArray();
	for(i=0; i<array.length/7; i++){
		index=i*7;
		x=array[index+3];
		y=array[index+4];
		z=array[index+5];
		used=0;
		if(usedArray.length>0){
			for(ii=0; ii<usedArray.length/3; ii++){
				indexUsed=ii*3;
				x3=usedArray[indexUsed];
				y3=usedArray[indexUsed+1];
				z3=usedArray[indexUsed+2];
				if(x==x3&&y==y3&&z==z3){
					used=1;
				}
			}
		}
		if(used==0){
			tempArray=newArray();
			tempDistanceArray=newArray();
			for(ii=0; ii<array.length/7; ii++){
				index2=ii*7;
				x2=array[index2+3];
				y2=array[index2+4];
				z2=array[index2+5];
				if(x==x2&&y==y2&&z==z2){			
					tempArray=Array.concat(tempArray,index2);
					tempDistanceArray=Array.concat(tempDistanceArray,array[index2+6]);
				}
			}
			if(tempArray.length>1){
				Array.getStatistics(tempDistanceArray,min,max,mean,std);
				for(ii=0; ii<tempDistanceArray.length; ii++){
					tempDist=tempDistanceArray[ii];
					index3=tempArray[ii];
					if(tempDist==min){
						cleanedUpConnections=Array.concat(cleanedUpConnections,array[index3],array[index3+1],array[index3+2],array[index3+3],array[index3+4],array[index3+5]);
						usedArray=Array.concat(usedArray,array[index3+3],array[index3+4],array[index3+5]);
					}
				}
			} else {
				cleanedUpConnections=Array.concat(cleanedUpConnections,array[index],array[index+1],array[index+2],array[index+3],array[index+4],array[index+5]);
				usedArray=Array.concat(usedArray,array[index+3],array[index+4],array[index+5]);
			}
		}
	}
	cleanedUpConnections2=newArray();
	for(i=0; i<cleanedUpConnections.length/6; i++){
		index=i*6;
		x=cleanedUpConnections[index];
		y=cleanedUpConnections[index+1];
		z=cleanedUpConnections[index+2];
		x1=cleanedUpConnections[index+3];
		y1=cleanedUpConnections[index+4];
		z1=cleanedUpConnections[index+5];
		distance=sqrt(pow(x*width-x1*width,2)+pow(y*height-y1*height,2)+pow(z*depth-z1*depth,2));
		if(distance!=0){
			cleanedUpConnections2=Array.concat(cleanedUpConnections2,cleanedUpConnections[index],cleanedUpConnections[index+1],cleanedUpConnections[index+2],cleanedUpConnections[index+3],cleanedUpConnections[index+4],cleanedUpConnections[index+5]);
		}
	}
	return cleanedUpConnections2;
}
/*
 * Create list of points that were not robustly connected
 */
function createNewPointList(robustConnectionsAnt,anteriorConnections){
	points=newArray();//This goes through and finds all objects that were not connected 
	row=0;
	for(i=0; i<anteriorConnections.length/9; i++){
		index=i*9;
		objectx=anteriorConnections[index];
		objecty=anteriorConnections[index+1];
		objectz=anteriorConnections[index+2];		
		if(robustConnectionsAnt[i]==0){
			points=Array.concat(points,objectx,objecty,objectz);
		}	
	}
	return points;
}


 
 /*
 * Find neighbors in the anterior direction
 */
function makeConnectionsShell(points,clusterSize,poles,scoring){
	totalConnectionsAnt=makeConnections(points,clusterSize,poles,1,scoring);
	totalConnectionsPost=makeConnections(points,clusterSize,poles,0,scoring);
	returnArray=Array.concat(totalConnectionsAnt,totalConnectionsPost);
	return returnArray;
}
 /*
 * Find neighbors in the anterior direction
 */
function makeConnectionsShellLoop(antpoints,postpoints,poles,scoring){
	totalConnectionsAnt=makeConnectionsLoop(antpoints,postpoints,postpoints.length/3,poles,1,scoring);
	totalConnectionsPost=makeConnectionsLoop(postpoints,antpoints,antpoints.length/3,poles,0,scoring);
	returnArray=Array.concat(totalConnectionsAnt,totalConnectionsPost);
	return returnArray;
}
/*
 * 
 */
function makeConnections(points,clusterSize,poles,antpost,scoring){
	totalConnections=newArray();
	for(i=0; i<100; i++){
		if(i==0){
			objects=points;
			neighbors=clusterNthObjects(15,objects);
			planeScores=planePenalty(objects,neighbors,poles,scoring);
			distanceScores=partnerDistancePenalty(objects,neighbors);
			poleScores=polePenalty(objects,neighbors,poles,antpost);
			connections=sortScores(15,planeScores,distanceScores,poleScores,neighbors,objects);
			reconnect=findPointsToReconnect(connections,objects);
			reconnectPartners=findPartnersToReconnect(connections,objects);
		} else{
			objects=reconnect;		
			neighbors=clusterNthObjectsLoop(reconnect,reconnectPartners);
			planeScores=planePenalty(objects,neighbors,poles,scoring);
			distanceScores=partnerDistancePenalty(objects,neighbors);
			poleScores=polePenalty(objects,neighbors,poles,antpost);
			connections=sortScores(reconnectPartners.length/3,planeScores,distanceScores,poleScores,neighbors,objects);
			reconnectPartners=findPartnersToReconnect(connections,objects);
			reconnect=findPointsToReconnect(connections,objects);
		}
		for(ii=0; ii<connections.length/6; ii++){
			index=ii*6; 
			index2=ii*3;
			if(connections[index]!=-1){
				totalConnections=Array.concat(totalConnections,objects[index2],objects[index2+1],objects[index2+2],connections[index],connections[index+1],connections[index+2],connections[index+3],connections[index+4],connections[index+5]);
			}
		}
		if(reconnect.length==0){
			i=100;
		}
	}
	return totalConnections;
}
/*
 * Loop version
 */
function makeConnectionsLoop(points,points2,clusterSize,poles,antpost,scoring){
	totalConnections=newArray();
	for(i=0; i<100; i++){
		if(i==0){
			objects=points;
			neighbors=clusterNthObjectsLoop(points,points2);
			planeScores=planePenalty(objects,neighbors,poles,scoring);
			distanceScores=partnerDistancePenalty(objects,neighbors);
			poleScores=polePenalty(objects,neighbors,poles,antpost);
			connections=sortScores(points2.length/3,planeScores,distanceScores,poleScores,neighbors,objects);
			reconnect=findPointsToReconnect(connections,objects);
			reconnectPartners=findPartnersToReconnect(connections,objects);
		} else{
			objects=reconnect;		
			neighbors=clusterNthObjectsLoop(reconnect,reconnectPartners);
			planeScores=planePenalty(objects,neighbors,poles,scoring);
			distanceScores=partnerDistancePenalty(objects,neighbors);
			poleScores=polePenalty(objects,neighbors,poles,antpost);
			connections=sortScores(reconnectPartners.length/3,planeScores,distanceScores,poleScores,neighbors,objects);
			reconnectPartners=findPartnersToReconnect(connections,objects);
			reconnect=findPointsToReconnect(connections,objects);
		}
		for(ii=0; ii<connections.length/6; ii++){
			index=ii*6; 
			index2=ii*3;
			if(connections[index]!=-1){
				totalConnections=Array.concat(totalConnections,objects[index2],objects[index2+1],objects[index2+2],connections[index],connections[index+1],connections[index+2],connections[index+3],connections[index+4],connections[index+5]);
			}
		}
		if(reconnect.length==0){
			i=100;
		}
	}
	return totalConnections;
}
/*
 * Find neighbors in the posterior direction
 */
/*
 *  Find connections that are robust in that the anterior and posterior connections are identical
 */
function keepRobustConnections(array1,array2){
	getVoxelSize(width,height,depth,unit);
	arrayReturn=newArray(array1.length/9);
	for(i=0; i<array1.length/9; i++){
		index=i*9;
		x=array1[index];//These are the objects coords.
		y=array1[index+1];
		z=array1[index+2];
		xant=array1[index+3];//These are the coords for what the object is connected to.
		yant=array1[index+4];
		zant=array1[index+5];
		distance=sqrt(pow(x*width-xant*width,2)+pow(y*height-yant*height,2)+pow(z*depth-zant*depth,2));
		for(ii=0; ii<array2.length/9; ii++){
			index2=ii*9;
			xantPost=array2[index2];
			yantPost=array2[index2+1];
			zantPost=array2[index2+2];
			if(xant==xantPost && yant==yantPost && zant==zantPost){
				xpost=array2[index2+3];
				ypost=array2[index2+4];
				zpost=array2[index2+5];
				if(x==xpost && y==ypost && z==zpost && distance<5){
					arrayReturn[i]=1;
				}else{
					arrayReturn[i]=0;
				}
			}
		}
		if(distance==0){
			arrayReturn[i]=0;
		}
	}
	return arrayReturn;
}
function crossProduct(array){
	vecA1=array[0]-(array[3]);
	vecA2=array[1]-(array[4]);
	vecA3=array[2]-(array[5]);
	vecB1=array[0]-(array[6]);
	vecB2=array[1]-(array[7]);
	vecB3=array[2]-(array[8]);
	vecC1=vecA2*vecB3-vecA3*vecB2;
	vecC2=vecA3*vecB1-vecA1*vecB3;
	vecC3=vecA1*vecB2-vecA2*vecB1;
	cp=newArray(vecC1,vecC2,vecC3);
	return cp;
}
function createPlaneEquation(array){
	xConstant=array[0];
	yConstant=array[1];
	zConstant=array[2];
	offset=(array[0]*(-1*array[3]))+(array[1]*(-1*array[4]))+(array[2]*(-1*array[5]));
	planeEquation=newArray(xConstant,yConstant,zConstant,offset);
	return planeEquation;
}
function distPointPlane(array){
	topTerm=abs((array[3]*array[0])+(array[4]*array[1])+(array[5]*array[2])+(array[6]));
	bottomTerm=sqrt(pow(array[3],2)+pow(array[4],2)+pow(array[5],2));
	distance=topTerm/bottomTerm;
	return distance;
}
function planePenalty(objects,neighbors,poles,scoring){
	getVoxelSize(width,height,depth,unit);
	clusterSize=neighbors.length/objects.length;
	arrayReturn=newArray();
	counter=0;
	for(i=0; i<objects.length/3; i++){
		index=i*3;
		point=newArray(objects[index],objects[index+1],objects[index+2]);
		array=Array.concat(poles,point);
		cp=crossProduct(array);
		array=Array.concat(cp,point);
		equation=createPlaneEquation(array);
		tempArrayReturn=newArray(clusterSize);
		index2=i*clusterSize*3;
		for(ii=0; ii<clusterSize; ii++){
			index3=index2+(ii*3);
			point2=newArray(neighbors[index3],neighbors[index3+1],neighbors[index3+2]);
			array=Array.concat(point2,equation);
			pointPlaneDistance=distPointPlane(array);
			if(scoring==2){
				arrayReturn=Array.concat(arrayReturn,pointPlaneDistance);
			}
			if(scoring==1){
				tempArrayReturn[ii]=pointPlaneDistance;
			}
			counter++;
		}
		if(scoring==1){
			rank=Array.rankPositions(tempArrayReturn);
			rank2=Array.rankPositions(rank);
			arrayReturn=Array.concat(arrayReturn,rank2);
		}
	}
	if(scoring==3){
		for(i=0; i<arrayReturn.length; i++){
			arrayReturn[i]=1;
		}
	}
	return arrayReturn;
}
function partnerDistancePenalty(objects,neighbors){
	getVoxelSize(width,height,depth,unit);
	clusterSize=neighbors.length/objects.length;
	arrayReturn=newArray((objects.length/3)*(clusterSize));
	counter=0;
	for(i=0; i<objects.length/3; i++){
		index=i*3;
		x1=objects[index]*width;//0,3,6,9,12...x coordinate for the object
		y1=objects[index+1]*height;//1,4,7,10,13...y coordinate for the object
		z1=objects[index+2]*depth;//2,5,8,11,14...z coordinate for the object
		index2=i*clusterSize*3;
		for(ii=0; ii<clusterSize; ii++){
			index3=index2+(ii*3);
			x2=neighbors[index3]*width;
			y2=neighbors[index3+1]*height;
			z2=neighbors[index3+2]*depth;
			distance=sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));
			arrayReturn[counter]=distance;
			counter++;
		}
	}
	return arrayReturn;
}
/*
 * Sort scores to find the lowest score. 
 * Return the components of the lowest score along with the xyz coordinates for the lowest score
 */
function sortScores(cluster,score1,score2,score3,neighbors,points){
	getVoxelSize(width,height,depth,unit);
	arrayReturn=newArray((score1.length/cluster)*6);//This array will return the XYZ and individual score components for every object along with a binary value indicating whether the object 
	scoreArray=newArray(score1.length);
	for(i=0; i<score1.length; i++){//This multiplies all the scores together to get a composite score
		scoreArray[i]=score1[i]*score2[i]*score3[i];
	}
	for(i=0; i<scoreArray.length/cluster; i++){//This loop goes through and connects everything based solely on its score
		index=i*cluster; 
		index2=i*6;
		index3=index*3;
		tempArray=newArray(cluster);
		for(ii=0; ii<cluster; ii++){//This fills a temporary arry with the scores for that objects
			tempArray[ii]=scoreArray[index+ii];
		}
		ranked=Array.rankPositions(tempArray);
		rankedScores=Array.rankPositions(ranked);
		for(ii=0; ii<cluster; ii++){
			if(rankedScores[ii]==1){
				arrayReturn[index2]=neighbors[((index+ii)*3)];
				arrayReturn[index2+1]=neighbors[((index+ii)*3)+1];
				arrayReturn[index2+2]=neighbors[((index+ii)*3)+2];
				arrayReturn[index2+3]=score1[index+ii];
				arrayReturn[index2+4]=score2[index+ii];
				arrayReturn[index2+5]=score3[index+ii];
			}
		}
	}
	deletePositions=newArray();
	for(i=0; i<arrayReturn.length/6; i++){//This loop goes back and finds anything that's received two connections
		index=i*6;
		index2=i*3;
		objectX=points[index2];
		objectY=points[index2+1];
		objectZ=points[index2+2];
		tempX=arrayReturn[index];//X coordinate of the object
		tempY=arrayReturn[index+1];//Y coordinate of the object
		tempZ=arrayReturn[index+2];//Z coordinate of the object
		tempScore=arrayReturn[index+4];//Distance between the object and its connected partner which is used to break ties
		duplicate=0;//This will stay 0 unless the search finds the same object more than once
		tempScoreSort=newArray();//An array to hold the scores of objects that are connected more than once
		tempPositionSort=newArray();//An array to hold the positions of the objects that are connected more than once. This should be the index
		for(ii=0; ii<arrayReturn.length/6; ii++){//This loop starts searching the partners for connections that were made more than once.
			index3=ii*6;
			tempX1=arrayReturn[index3];//X coordinate of the connected object
			tempY1=arrayReturn[index3+1];//Y coordinate of the connected object
			tempZ1=arrayReturn[index3+2];//Z coordinate of the connected object
			tempScore1=arrayReturn[index3+4];//Distance between the connected object and its partner
			if(ii!=i){//This prevents the search from finding the connected object itself
				if(tempX==tempX1 && tempY==tempY1 && tempZ==tempZ1){//This is only satisfied if it found the same object
					duplicate=1;//This lets the loop know that it found the same connected object
					tempScoreSort=Array.concat(tempScoreSort,tempScore1);//This adds the score to the score array
					tempPositionSort=Array.concat(tempPositionSort,index3);//This adds the position of multiple connected objects to the position array
				}
			}
		}
		if(duplicate==1){//If a multiple connected object was found than add the initial object to the array
			tempScoreSort=Array.concat(tempScoreSort,tempScore);
			tempPositionSort=Array.concat(tempPositionSort,index);
			rank=Array.rankPositions(tempScoreSort);
			rankedScores=Array.rankPositions(rank);
			for(ii=0; ii<rankedScores.length; ii++){
				if(rankedScores[ii]!=0){
					deleteIndex=tempPositionSort[ii];
					deletePositions=Array.concat(deletePositions,deleteIndex);
				}
			}
		}
		duplicate=0;
	}
	for(i=0; i<deletePositions.length; i++){
		index=deletePositions[i];
		arrayReturn[index]=-1;
		arrayReturn[index+1]=-1;
		arrayReturn[index+2]=-1;
		arrayReturn[index+3]=-1;
		arrayReturn[index+4]=-1;
		arrayReturn[index+5]=-1;
	}
	return arrayReturn;
}
/*
 * Find n closest objects to object of interest.
 * Cluster is the number of objects to include in the cluster.
 * array is an array of coordinates in xyz format.
 * The function returns an array of xyz coordinates that contain the neighbors for
 * every object in array. For the array to be deconvolved the array needs to be paired 
 * with the same cluster size used to create the array.
 */ 
function clusterNthObjects(cluster,array){
	getVoxelSize(width,height,depth,unit);
	clusterCoordinatesArray=newArray();//This will be the finall array that has all the coordinates for the clustered objects
	for(i=0; i<array.length/3; i++){//This marches through all the objects of interest
		index=i*3;//0,3,6,9,12...This is the index for pulling the coordinates for the object from the coordinatesArray
		x1=array[index]*width;//0,3,6,9,12...x coordinate for the object
		y1=array[index+1]*height;//1,4,7,10,13...y coordinate for the object
		z1=array[index+2]*depth;//2,5,8,11,14...z coordinate for the object
		distanceArray=newArray(array.length/3);//This array contains all the distances relative to the object of interest
		for(ii=0; ii<array.length/3; ii++){//This marches through all objects a second time
			index2=ii*3;//0,3,6,9,12...This index is the index for pulling the objects the coordinates for the object from the coordinatesArray
			x2=array[index2]*width;//0,3,6,9,12...x coordinate for the 2nd object
			y2=array[index2+1]*height;//1,4,7,10,13...y coordinate for the 2nd object
			z2=array[index2+2]*depth;//2,5,8,11,14...z coordinate for the 2nd object
			distance=sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));//This is the distance between the first and second object
			distanceArray[ii]=distance;//This assigns the distance between the two objects to the distanceArray
			
		}
		ranked=Array.rankPositions(distanceArray);//This ranks the distanceArray so that a new array is created that contains the relative positions of the 2nd objects relative to the first object
		rankedDistanceArray=Array.rankPositions(ranked);
		tempClusterCoordinatesArray=newArray(cluster*3);//This creates a new array that has enough space to hold the xyz coordinates for nth closest objects
		counter=0;
		for(ii=0; ii<array.length/3; ii++){//This loop marches through the rankedDistanceArray looking for positions that are less than cluster size
			value=rankedDistanceArray[ii];
			if(value<cluster){//If the rankedDistance is less than the cluster size then extract the spatial coordinates
				index2=counter*3;
				tempClusterCoordinatesArray[index2]=array[(ii*3)];
				tempClusterCoordinatesArray[index2+1]=array[(ii*3)+1];
				tempClusterCoordinatesArray[index2+2]=array[(ii*3)+2];
				counter++;
			}
		}
		clusterCoordinatesArray=Array.concat(clusterCoordinatesArray,tempClusterCoordinatesArray);
	}
	return clusterCoordinatesArray;
}
/*
 * Find n closest objects to object of interest.
 * Cluster is the number of objects to include in the cluster.
 * array is an array of coordinates in xyz format.
 * The function returns an array of xyz coordinates that contain the neighbors for
 * every object in array. For the array to be deconvolved the array needs to be paired 
 * with the same cluster size used to create the array.
 */ 
function clusterNthObjectsLoop(array,array2){
	getVoxelSize(width,height,depth,unit);
	cluster=array2.length/3;
	clusterCoordinatesArray=newArray();//This will be the finall array that has all the coordinates for the clustered objects
	for(i=0; i<array.length/3; i++){//This marches through all the objects of interest
		index=i*3;//0,3,6,9,12...This is the index for pulling the coordinates for the object from the coordinatesArray
		x1=array[index]*width;//0,3,6,9,12...x coordinate for the object
		y1=array[index+1]*height;//1,4,7,10,13...y coordinate for the object
		z1=array[index+2]*depth;//2,5,8,11,14...z coordinate for the object
		distanceArray=newArray(array2.length/3);//This array contains all the distances relative to the object of interest
		for(ii=0; ii<array2.length/3; ii++){//This marches through all objects a second time
			index2=ii*3;//0,3,6,9,12...This index is the index for pulling the objects the coordinates for the object from the coordinatesArray
			x2=array2[index2]*width;//0,3,6,9,12...x coordinate for the 2nd object
			y2=array2[index2+1]*height;//1,4,7,10,13...y coordinate for the 2nd object
			z2=array2[index2+2]*depth;//2,5,8,11,14...z coordinate for the 2nd object
			distance=sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));//This is the distance between the first and second object
			distanceArray[ii]=distance;//This assigns the distance between the two objects to the distanceArray
		}
		ranked=Array.rankPositions(distanceArray);//This ranks the distanceArray so that a new array is created that contains the relative positions of the 2nd objects relative to the first object
		rankedDistanceArray=Array.rankPositions(ranked);
		tempClusterCoordinatesArray=newArray(cluster*3);//This creates a new array that has enough space to hold the xyz coordinates for nth closest objects
		counter=0;
		for(ii=0; ii<array2.length/3; ii++){//This loop marches through the rankedDistanceArray looking for positions that are less than cluster size
			value=rankedDistanceArray[ii];
			if(value<cluster){//If the rankedDistance is less than the cluster size then extract the spatial coordinates
				index2=counter*3;
				tempClusterCoordinatesArray[index2]=array2[(ii*3)];
				tempClusterCoordinatesArray[index2+1]=array2[(ii*3)+1];
				tempClusterCoordinatesArray[index2+2]=array2[(ii*3)+2];
				counter++;
			}
		}
		clusterCoordinatesArray=Array.concat(clusterCoordinatesArray,tempClusterCoordinatesArray);
	}
	return clusterCoordinatesArray;
}
function findPointsToReconnect(connected,points){
	getVoxelSize(width,height,depth,unit);
	reconnectArray=newArray();
	for(i=0; i<connected.length/6; i++){//This loop goes through and connects everything based solely on its score
		index=i*6;
		index2=i*3;
		objectx=points[index2];
		objecty=points[index2+1];
		objectz=points[index2+2];		
		if(connected[index]==-1){
			reconnectArray=Array.concat(reconnectArray,objectx,objecty,objectz);
		}
	}
	return reconnectArray;
}
/*
 * Find possible partners to reconnect to
 */
function findPartnersToReconnect(connected,points){
	getVoxelSize(width,height,depth,unit);
	possibleArray=newArray();
	for(i=0; i<connected.length/6; i++){//This loop goes through and connects everything based solely on its score
		index=i*6;
		index2=i*3;
		objectx=points[index2];
		objecty=points[index2+1];
		objectz=points[index2+2];
		possible=0;
		for(ii=0; ii<connected.length/6; ii++){
			index3=ii*6;
			queryx=connected[index3];
			queryy=connected[index3+1];
			queryz=connected[index3+2];
			if(objectx==queryx && objecty==queryy && objectz==queryz){
				possible=1;
			}
		}
		if(possible==0){
			possibleArray=Array.concat(possibleArray,objectx,objecty,objectz);
		}
	}
	return possibleArray;
}
/*
 * 
 */
function interactiveCursor(boxWidth,boxDepth,currentImage,binaryImage){
	if (getVersion>="1.37r"){
		setOption("DisablePopupMenu", true);
	}
	getPixelSize(unit, pixWidth, pixHeight);
	shift=1;
	ctrl=2; 
	rightButton=4;
	alt=8;
	leftButton=16;
	insideROI = 32;
	x2=-1; y2=-1; z2=-1; flags2=-1;
	stop=false;
	bbArray=newArray();
	while (!stop) {
		selectImage(currentImage);
		resetMinAndMax();
		getCursorLoc(x, y, z, flags);
		overlay=false;
		if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) {
			s = " ";
			if (flags&leftButton!=0){
				hyperStack=Stack.isHyperstack;
				if(hyperStack==true){
					Stack.getPosition(channel, slice, frame);
				}else{
					slice=z;
				}
				selectImage(binaryImage);
				array=centerBox(x,y,slice,boxWidth,boxDepth,1);
				centerx=parseFloat(array[0]);
				centery=parseFloat(array[1]);
				centerz=parseFloat(array[2]);
				upperXFinal=parseFloat(array[3]);
				upperYFinal=parseFloat(array[4]);
				boxWidthPix=parseFloat(array[5]);
				boxDepthPix=parseFloat(array[6]);
				overlay=true;
				bbArray=Array.concat(bbArray,centerx,centery,centerz);
				wait(250);
			}
			if (flags&rightButton!=0){
				stop=true;
			}
			if (flags&shift!=0){
				waitForUser("shit");
			}
	
		}
		x2=x; y2=y; z2=z; flags2=flags;
		if(bbArray.length>1 && overlay==true){
			selectImage(currentImage);
			if(hyperStack==true){
				Stack.setSlice(centerz);
			}else{
				setSlice(centerz);
			}
			setColor("green");
			Overlay.drawRect(upperXFinal,upperYFinal,boxWidthPix,boxWidthPix);
			if(hyperStack==true){
				Overlay.setPosition(1,centerz,1);
			}else{
				Overlay.setPosition(centerz);
			}
			Overlay.setPosition(centerz);
			Overlay.show();
			wait(100);
		} 
		wait(100);
	}
	if (getVersion>="1.37r"){
		setOption("DisablePopupMenu", false);
	}
	return bbArray;
}
/*
 * Center the box on the highest intensity voxel 
 */
function centerBox(x,y,slice,boxWidth,boxDepth,zSearch){
	title2=getImageID();
	getVoxelSize(width,height,depth,unit);
	totalSlices=nSlices;
	boxWidthPix=round(boxWidth/width);
	boxDepthPix=round(boxDepth/depth);
	if(boxWidthPix%2==0){
		boxWidthPix=boxWidthPix+1;
	}
	if(boxDepthPix%2==0){
		boxDepthPix=boxDepthPix+1;
	}
	setSlice(slice);
	makeRectangle(x-((boxWidthPix-1)/2),y-((boxWidthPix-1)/2),boxWidthPix,boxWidthPix);
	getSelectionBounds(upperX, upperY, void, void);
	//waitForUser("center");
	if(zSearch==1){
		startSlice=slice-((boxDepthPix-1)/2);
		startRemainder=0;
		stopRemainder=0;
		frontTruncated=false;	
		if(startSlice<=0){
			startRemainder=abs(startSlice-1);
			startSlice=1;
			frontTruncated=true;
		}
		stopSlice=slice+((boxDepthPix-1)/2);
		backTruncated=false;
		if(stopSlice>totalSlices){
			stopRemainder=abs(stopSlice-totalSlices);
			stopSlice=totalSlices;
			backTruncated=true;
		}
		rangeString=toString(startSlice)+"-"+toString(stopSlice);
		run("Duplicate...", "title=[Temp] duplicate range=&rangeString");
		bits=bitDepth;
		temp=getImageID();
		radiusX=getWidth;
		radiusY=getHeight;
		radiusZ=nSlices;
		if(frontTruncated==true){
			if(bits==8){
				newImage("TempAdd", "8-bit black", radiusX, radiusY, startRemainder);
			}
			if(bits==16){
				newImage("TempAdd", "16-bit black", radiusX, radiusY, startRemainder);
			}
			run("Concatenate", "  title=[Concatenated Stacks] image_1=TempAdd image_2=Temp");
			rename("Temp");
			temp=getImageID();
		}
		if(backTruncated==true){
			if(bits==8){
				newImage("TempAdd", "8-bit black", radiusX, radiusY, stopRemainder);
			}
			if(bits==16){
				newImage("TempAdd", "16-bit black", radiusX, radiusY, stopRemainder);
			}
			run("Concatenate", "  title=[Concatenated Stacks] image_1=Temp image_2=TempAdd");
			rename("Temp");
			temp=getImageID();
		}
		run("3D Fast Filters","filter=MaximumLocal radius_x_pix=&radiusX radius_y_pix=&radiusX radius_z_pix=&radiusZ Nb_cpus=4");
		temp2=getImageID();
		Stack.getStatistics(voxels,mean,min,max,std);
		coordsArray=newArray();
		midX=((radiusX-1)/2);
		midY=((radiusY-1)/2);
		midZ=((radiusZ-1)/2)+1;
		for(i=0; i<radiusZ; i++){
			setSlice(i+1);
			for(ii=0; ii<radiusX; ii++){
				for(iii=0; iii<radiusX; iii++){
					intensity=getPixel(ii,iii);
					if(intensity==max){
						coordsArray=Array.concat(coordsArray,ii,iii,i+1);
					}
				}
			}
		}
		if(coordsArray.length>3){
			distanceArray=newArray();
			for(i=0; i<coordsArray.length/3; i++){
				x1=coordsArray[(i*3)+0];
				y1=coordsArray[(i*3)+1];
				z1=coordsArray[(i*3)+2];
				distance=sqrt(pow(x1-midX,2)+pow(y1-midY,2)+pow(z1-midZ,2));
				distanceArray=Array.concat(distanceArray,distance);
			}
			Array.getStatistics(distanceArray,min,max,mean,std);
			for(i=0; i<distanceArray.length; i++){
				if(min==distanceArray[i]){
					centerX=upperX+coordsArray[(i*3)+0];
					centerY=upperY+coordsArray[(i*3)+1];
					if(coordsArray[(i*3)+2]==midZ){
						centerZ=slice;
					}
					if(coordsArray[(i*3)+2]<midZ){
						centerZ=slice-(midZ-coordsArray[(i*3)+2]);
					}
					if(coordsArray[(i*3)+2]>midZ){
						centerZ=slice+(coordsArray[(i*3)+2]-midZ);
					}
				}
			}
	
		}else{
			centerX=upperX+coordsArray[0];
			centerY=upperY+coordsArray[1];	
			if(coordsArray[2]==midZ){
				centerZ=slice;
			}
			if(coordsArray[2]<midZ){
				centerZ=slice-(midZ-coordsArray[2]);
			}
			if(coordsArray[2]>midZ){
				centerZ=slice+(coordsArray[2]-midZ);
			}		
		}
		selectImage(temp);
		close();
		selectImage(temp2);
		close();		
	}
	if(zSearch==0){
		run("Duplicate...", "title=[Temp]");
		//waitForUser("Duplicate");
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		temp=getImageID();
		radiusX=getWidth;
		radiusY=getHeight;
		intensityArray=newArray();
		coordsArray=newArray();
		midX=((radiusX-1)/2);
		midY=((radiusY-1)/2);
		for(i=0; i<radiusX; i++){
			for(ii=0; ii<radiusY; ii++){
				intensity=getPixel(i,ii);
				if(intensity==max){
					coordsArray=Array.concat(coordsArray,i,ii,1);
				}
			}
		}
		if(coordsArray.length>3){
			distanceArray=newArray();
			for(i=0; i<coordsArray.length/3; i++){
				x1=coordsArray[(i*3)+0];
				y1=coordsArray[(i*3)+1];
				distance=sqrt(pow(x1-midX,2)+pow(y1-midY,2));
				distanceArray=Array.concat(distanceArray,distance);
			}
			Array.getStatistics(distanceArray,min,max,mean,std);
			for(i=0; i<distanceArray.length; i++){
				if(min==distanceArray[i]){
					centerX=upperX+coordsArray[(i*3)+0];
					centerY=upperY+coordsArray[(i*3)+1];
				}
			}
			
		}else{
			centerX=upperX+coordsArray[0];
			centerY=upperY+coordsArray[1];
		}
		selectImage(temp);
		close();
		centerZ=slice;		
	}
	upperXFinal=centerX-((boxWidthPix-1)/2);
	upperYFinal=centerY-((boxWidthPix-1)/2);
	selectImage(title2);
	run("Select None");
	returnArray=newArray(centerX,centerY,centerZ,upperXFinal,upperYFinal,boxWidthPix,boxDepthPix);
	return returnArray;
}
/*
 * 
 */
function modifyBinaryImage(image,modification,coords,boxSize,selection){
	selectImage(image);
	for(i=0; i<coords.length/3; i++){
		setSlice(coords[(i*3)+2]);
		if(selection=="square"){
			makeRectangle(coords[(i*3)+0],coords[(i*3)+1],boxSize,boxSize);
		}
		if(modification=="delete"){
			run("Clear","slice");
		}else if(modification=="add"){
			if(bitDepth==8){
				setColor(255);
				fill();
			}
			if(bitDepth==16){
				setColor(65535);
				fill();
			}
			run("Select None");
		}
	}
}
/*
 * 
 */
function parseArray(input){
	output=newArray();
	for(i=0; i<input.length; i++){
		temp=parseFloat(input[i]);
		output=Array.concat(output,temp);
	}
	return output;
}
