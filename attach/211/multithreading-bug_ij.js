// demonstrate possible IJ multithreading bug
// multithreaded image analysis in javascript
// Thomas Siegmund and Sunaina Prasad @ Evotec AG

var currentTime = new Date();
var startTime = currentTime.getTime();

importPackage(java.lang);
importPackage(java.io);
importClass(java.io.File);
importClass(java.io.FileWriter);
importClass(java.io.BufferedWriter);
importClass(Packages.java.util.concurrent.atomic.AtomicInteger);

importClass(Packages.ij.IJ);
importClass(java.awt.Color);
importClass(java.awt.Font); 
importClass(Packages.ij.gui.Overlay);
importClass(Packages.loci.plugins.BF);
importClass(Packages['loci.plugins.in.ImporterOptions']);

// new stuff needed since 2014/06
importPackage(Packages.ij.measure);
importPackage(Packages.ij.plugin);
importPackage(Packages.ij.plugin.filter);
importPackage(Packages.ij.process);
importPackage(Packages.ij.io);
importClass(Packages.ij.ImagePlus);
importPackage(Packages.ij.plugin.frame);
importPackage(Packages.ij.gui);

var phdnThresholder = "Percentile dark";
var trtThresholder = "Triangle dark";
var nucleiThresholder = "Otsu dark";
var nucleiNoiseTolerance = 10;
var nucleiSigma = 1.2;
var nucleusMinSize = 20;
var nucleusMaxSize = 400;
var nucleusMinRoundness = 0.3;
var nucleiThresholder = "Otsu dark";
var cytoplasmMinSize = 400;
var cytoplasmMaxSize = 500000;
var clusterMinRoundness = 0.01; //original 0.1

var	firstImgProcessed = false;
var	lowerDThreshold=0;
var	upperDThreshold=0;
var	lowerCytoThreshold =0;
var upperCytoThreshold=0;

// restrict analyis to central oval
var restrictSelection = false;
var boarderWidth = 400;
var testing = false;

// remove extra bright nuclei which we think are primarily debris
var removeBrightNuclei = true;
var brightNucleiThresholder = "Minimum dark";

var dc = DirectoryChooser("Select an input folder...");
var outformat = "tif";
var path = dc.getDirectory();
path = path.slice(0,path.length()-1);
Packages.ij.macro.Interpreter.batchMode = true; 

multithreading = true;

// per image worker function
function processFile(i) {

	var image = path+"/"+imagefiles[i];
	IJ.log("opening: "+image);

	IJ.run("Clear Results");
	// load image, split channels
	var options = new ImporterOptions();
	options.setId(image);
	// do not use autoscale as we measure intensities
	options.setAutoscale(false);
	options.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE);
	var img = BF.openImagePlus(options)[0];
	var ip = img.getProcessor();
	ip.setColor(Color.red);
	var nchannels = img.getNChannels();
	var re = /[^\.]*/;
	var currentImageName = re.exec(imagefiles[i]);
	var channels = ChannelSplitter.split(img);
//	var stack = img.getStack();
	
	var dapiImg = channels[0];
	var dapiIp = dapiImg.getProcessor();
	IJ.run(dapiImg, "Grays", "");
	IJ.run(dapiImg, "Subtract...", "value=" + lowvalues[0]);
	IJ.run(dapiImg, "Multiply...", "value=" + 65536 / (highvalues[0] - lowvalues[0]));
	IJ.resetMinAndMax(dapiImg);
	dapiImg.setTitle("Dapi");
	if(testing) {
		dapiImg.show();
	}

	var phdnImg = channels[1];
	var phdnIp = phdnImg.getProcessor();
	IJ.run(phdnImg, "Grays", "");
	IJ.run(phdnImg, "Subtract...", "value=" + lowvalues[2]);
	IJ.run(phdnImg, "Multiply...", "value=" + 65536 / (highvalues[2] - lowvalues[2]));
	IJ.resetMinAndMax(phdnImg);
	phdnImg.setTitle("Phdn");
	
	if(testing) {
		phdnImg.show();
	}

  var trtImg = channels[2];
	var trtIp = trtImg.getProcessor();
	IJ.run(trtImg, "Grays", "");
	IJ.run(trtImg, "Subtract...", "value=" + lowvalues[1]);
	IJ.run(trtImg, "Multiply...", "value=" + 65536 / (highvalues[1] - lowvalues[1]));
	IJ.resetMinAndMax(trtImg);
	trtImg.setTitle("trt");
	
	if(testing) {
		trtImg.show();
	}

	if(channels.length == 4) {
	   	var smdImg = channels[3];
		var smdIp = smdImg.getProcessor();
		IJ.run(smdImg, "Grays", "");
		IJ.run(smdImg, "Subtract...", "value=" + lowvalues[1]);
		IJ.run(smdImg, "Multiply...", "value=" + 65536 / (highvalues[1] - lowvalues[1]));
		IJ.resetMinAndMax(smdImg);
		smdImg.setTitle("SMD");
		if(testing) {
			smdImg.show();
		}
	}

  img.changes = false;
  img.close();
  img = null;
	 	
	// analyze total nuclei
	var rt_NucleiCount = new ResultsTable();
	rt_NucleiCount.reset();

	IJ.setThreshold(dapiImg, 18461, 61697);

	IJ.run(dapiImg, "Convert to Mask", "");	
	// here we get rid of most of the circles left over from the over-bright nuclei
	IJ.run(dapiImg, "Erode", "");

    var particle_analyzer = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW, 64, rt_NucleiCount, nucleusMinSize, nucleusMaxSize, nucleusMinRoundness, 1);
	particle_analyzer.analyze(dapiImg);
	var NucleiCount= rt_NucleiCount.getCounter();
	
	// use palloidin image to create cell mask 
	IJ.run(phdnImg, "Gaussian Blur...", "sigma=1");

  IJ.setThreshold(phdnImg, 6544, 65535);

	var phdnRoi = ThresholdToSelection.run(phdnImg);
	
	var cytoplasmImg = new Duplicator().run(trtImg);
	cytoplasmImg.setTitle("cytoplasm mask");
	if(testing) {
		cytoplasmImg.show();
	}
	phdnImg.setRoi(phdnRoi);
	trtImg.setRoi(phdnRoi);
	cytoplasmImg.setRoi(phdnRoi);

	if (restrictSelection== true){
	IJ.run(cytoplasmImg, "Clear Outside", "");
	}
	
	if(channels.length == 4) {
		smdImg.setRoi(phdnRoi);	
	}

	IJ.setThreshold(cytoplasmImg, 771, 65535);
	IJ.run(cytoplasmImg, "Convert to Mask", "");
	var trtRoi = ThresholdToSelection.run(cytoplasmImg);

	var fibroblastSmadIntensity = -1;
	
	// measure SMD in fibroblast nuclei
	fibroBlastNucleiMeanSmadIIntensity = -1;

	var rt_cytoplasmObjectCount = new ResultsTable();
	rt_cytoplasmObjectCount.reset();
	var manager = new RoiManager(true); // set to true to hide manager
  var cytoplasmPA = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW + ParticleAnalyzer.ADD_TO_MANAGER, 64, rt_cytoplasmObjectCount, cytoplasmMinSize, cytoplasmMaxSize, clusterMinRoundness, 1);
	cytoplasmPA.setRoiManager(manager);
	cytoplasmPA.analyze(cytoplasmImg);
	var rois = manager.getRoisAsArray();

	// using new array to store rois and closing the manager
	var roi_array= new Array();
	for (var i = 0; i < rois.length; i++) {
		roi_array[i]= rois[i];		
	}
	manager.runCommand("Delete");
	manager.close();

  IJ.setThreshold(cytoplasmImg, 1, 65535); // 
	var cytoplasmRoi = ThresholdToSelection.run(cytoplasmImg);

  dapiImg.deleteRoi();

	var rt_cmcsSMD = new ResultsTable();

  // loop over all cytoplasm objects
	for (var i = 0; i < roi_array.length; i++) { 
		roi_array[i].setStrokeColor(red); 
		trtImg.setRoi(roi_array[i]);
		dapiImg.setRoi(roi_array[i]);
		trtImg.deleteRoi();

    // now count and measure nuclei within
		var nucleiMeasurements = Measurements.AREA + Measurements.CIRCULARITY;
		var rt_nuclei = new ResultsTable();
		rt_nuclei.reset();

    // this particle analyzer stuff is strange. Output in the results table depends on "SHOW" settings unexpectedly. No output if set to "SHOW_NONE"
		var nucleiPA = new ParticleAnalyzer(ParticleAnalyzer.SHOW_RESULTS + ParticleAnalyzer.ADD_TO_MANAGER + ParticleAnalyzer.CLEAR_WORKSHEET, nucleiMeasurements, rt_nuclei, nucleusMinSize, nucleusMaxSize, nucleusMinRoundness, 1);
	
		var nuclei_manager = new RoiManager(true); // set to true to hide manager
		nucleiPA.setRoiManager(nuclei_manager);
	 	nucleiPA.analyze(dapiImg);
		var nNuclei = rt_nuclei.getCounter();
		
		// this is to obtain the nuclei within the cytoplasm
    // remove this line and the error seems to be gone
		var nucleiRois = nuclei_manager.getRoisAsArray();

    if(nNuclei != nuclei_manager.getCount()) {
			// this happens only if the script runs multithreaded. ???
			IJ.log(" ERROR: number of nuclei != number of ROIS !!!!: "+image+" "+i+", nNuclei "+nNuclei+", nucleiRoisLength "+nucleiRois.length);
		}
		
		// more image analysis on the nucleiRois to happen here
		
	}
	//nuclei_manager.runCommand("Delete");	
	nuclei_manager.close();

	cytoplasmImg.changes = false;
	cytoplasmImg.close();
	cytoplasmImg = null;

  System.gc();
	return(0);
	
} // end of processFile


// start of "main"
var red = new Color(1, 0, 0);
var orange = new Color(1, 0.43, 0);

IJ.run("Clear Results");
IJ.run("Close All");
IJ.run("Input/Output...", "jpeg=100 gif=-1 file=.txt use copy_row save_column save_row");

var dir = new File(path);
var imagefiles = new Array();
var files = dir.list();
if(files.length > 0) {

	var lowvalues = new Array();
	var highvalues = new Array();
	lowvalues[0] = 59.22;
	highvalues[0] = 1814.61;
	lowvalues[1] = 195.48;
	highvalues[1] = 1369.35;	
	lowvalues[2] = 42.039;
	highvalues[2] = 2490.06;

	// build file list
	for (var j = 0; j < files.length; j++) {
		var re = /.+(jpg|JPG|png|PNG|TIFF?|tiff?)$/;
		if(re.test(files[j])) {
			imagefiles.push(files[j]);
		}
	}

	for (var j = 0; j < 3; j++) {
	imagefiles.push('test_image.tiff');
	}

	processFile(0);
	
	if(multithreading) {
		multithreader(processFile, 1, imagefiles.length - 1);
	} else {
		for(i = 1; i < imagefiles.length; i++) {
		processFile(i);
		}
	}
}

Packages.ij.macro.Interpreter.batchMode = false; 
var currentTime = new Date();
var endTime = currentTime.getTime();
IJ.log("time: "+((endTime-startTime)/1000)+" seconds");
// end of "main"

// multithreading framework from
// http://fiji.sc/Javascript_Scripting#Multithreaded_Image_Processing_in_Javascript
function multithreader(fun, start, end) {
	var threads = new java.lang.reflect.Array.newInstance(java.lang.Thread, Runtime.getRuntime().availableProcessors());
	IJ.log("threads: "+Runtime.getRuntime().availableProcessors());
	var ai = new AtomicInteger(start);
	// Prepare arguments: all other arguments passed to this function
	// beyond the mandatory arguments fun, start and end:
	var args = new Array();
	var b = 0;
	for (var a = 3; a < arguments.length; a++) {
		args[b] = arguments[a];
		IJ.log("  argument " + (b+1) + " is " + args[b]);
		b++;
	}
	var body = {
		run: function() {
			for (var i = ai.getAndIncrement(); i <= end; i = ai.getAndIncrement()) {
				// Execute the function given as argument,
				// passing to it all optional arguments:
				fun(i, args);
			}
		}
	}
	// start all threads
	for (var i = 0; i < threads.length; i++) {
		threads[i] = new Thread(new Runnable(body)); // automatically as Runnable
		threads[i].start();
	}
	// wait until all threads finish
	for (var i = 0; i < threads.length; i++) {
		threads[i].join();
	}
}

