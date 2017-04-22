list=getFileList("D:\\cell mask\\control\\fused\\");
for(i=0;i<list.length;i++)
{
open("D:\\cell mask\\control\\fused\\" + list[i]);
imgName=getTitle(); 
//run("Threshold...");
setThreshold(70, 255);
setOption("BlackBackground", false);
run("Make Binary", "thresholded remaining black");
run("Smooth");
run("Invert");
run("Measure");
run("Close");
}

