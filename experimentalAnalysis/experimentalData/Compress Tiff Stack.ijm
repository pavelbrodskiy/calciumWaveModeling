input = "F:\\calcium waves simulation\\rawDiscData2\\Raw\\";
output = "F:\\calcium waves simulation\\rawDiscData2\\Compressed\\";

setBatchMode(true);

list = getFileList(input);

for (i = 0; i < list.length; i++) {


open(input + list[i]);
run("Enhance Contrast", "saturated=0.35");
run("8-bit");
newName = replace(list[i], ".tif", ".avi");
run("AVI... ", "frame=30 save=[F:\\calcium waves simulation\\rawDiscData2\\Compressed\\"+newName+"]");

}

setBatchMode(false);