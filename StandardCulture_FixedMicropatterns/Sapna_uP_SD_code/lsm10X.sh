#!/bin/bash
# running ilastik segmentation on 96wellplate - each well a separate folder
#trainproject: path corresponding to ilastik project that was used to train a sample dataset
#inputfiles: path corresponding to subfolders with sample images(*_ch3.tif) that are to be segmented



trainproject=/Volumes/Li/FigureData/20200302smad/allData/ch3.ilp

inputfiles=(
/Volumes/Li/FigureData/20200302smad/allData/processedData/*)

for f in ${inputfiles[@]}
do
f1=( 
$f/*_ch3.tif 
)

for ii in ${f1[@]}
do
	/Users/suprna/Downloads/ilastik-1.3.2post1-OSX.app/Contents/ilastik-release/run_ilastik.sh --headless --project=$trainproject $ii  --export_source="Simple Segmentation"
done	
done


