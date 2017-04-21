//Olexandra Ovsiy
//Macrofile for batch super-resolution analysis of two-color images
//using QuickPALM plugin
//script runs on linux without GUI display
//Linux command line: ./fiji --headless -batch test_cluster.ijm

setBatchMode(true);

dir = '/ifs/home/oo370/fiji/myTest/';
print(dir);
list = getFileList(dir);


Array.print(list);
for (i=0; i<list.length; i++) {
open(dir+list[i]);
run("Split Channels");
selectWindow("C2-"+list[i]);
saveAs("Tiff", dir+"C2-"+list[i]);
run("Analyse Particles", "minimum=2 maximum=4 image=127 smart online stream file=[/ifs/home/oo370/fiji/myTest/Particles Table.xls] pixel=20 accumulate=0 update=100 _image=imgNNNNNNNNN.tif start=0 in=50 _minimum=0 local=20 _maximum=1000 threads=50");
selectWindow("C2-"+(i+1)+" Reconstruction");
saveAs("Tiff", dir +"C2-"+(i+1)+" Reconstruction");
close();
}