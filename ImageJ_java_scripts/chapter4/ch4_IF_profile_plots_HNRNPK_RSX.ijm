//This script used for HNRNPK Figure 4-24 and associated Supplementary Figure S4-7 (Appendix C)
//Also used for TRIP12 Figure 4-25 and associated Supplementary Figure S4-8 (Appendix C)
//Raw data is .czi images and metadata from Zeiss confocal LSM780 and Zeiss ZEN software
//images are in single Z-plane 
//open original 16 bit image as Composite. Can alter contrast to view, but does not change pixel values. 
//For single line profile. All channel profiles saved in single .txt file
//manual input required: open image + adjust contrast + select line location + select crop
original=getTitle();
Stack.setActiveChannels("0101");
setTool("line");
waitForUser("set line");
roiManager("Add");
run("Split Channels");
selectWindow("C2-" + original);
rename("rsx");
run("Yellow");
roiManager("Select", 0);
run("Plot Profile");
Plot.setStyle(0, "black,none,1.0,Line");
//resetMinAndMax();
selectWindow("rsx");
waitForUser("adjust contrast");
selectWindow("C4-" + original);
rename("protein");
roiManager("Select", 0);
run("Plot Profile");
Plot.setStyle(0, "green,none,1.0,Line");
//resetMinAndMax();
selectWindow("protein");
waitForUser("adjust contrast");
selectWindow("C1-" + original);
close();
//rename("fibrillarin");
//run("Magenta");
//resetMinAndMax();
//waitForUser("adjust contrast");
selectWindow("C3-" + original);
close();
//rename("dapi");
//resetMinAndMax();
//waitForUser("adjust contrast");
run("Merge Channels...", "c2=protein c7=rsx create");
run("RGB Color");
selectWindow("Composite");
close();
rename("forline");
roiManager("Select", 0);
run("Flatten");
setTool("rectangle");
waitForUser("select for crop");
run("Crop");
selectWindow("forline");
close();
//rename(original + "_lineimage.tif");
//saveAs("tiff", "/Users/.../"+getTitle());
//close();
selectWindow("Plot of protein");
Plot.addFromPlot("Plot of rsx", 0);
Plot.setStyle(1, "black,none,1.0,Line");
selectWindow("Plot of rsx");
close();
//Plot.showValues();
selectWindow("Plot of protein");
rename(original+"_profile");
run("Images to Stack", "name=Stack title=[] use");
Stack.swap(1,2);
run("Make Montage...", "columns=2 rows=1 scale=1");
rename(original + "_profile.tif");
waitForUser("check");
saveAs("tiff", "/Users/.../"+getTitle());
roiManager("Select", 0);
roiManager("Delete");
run("Close All");

