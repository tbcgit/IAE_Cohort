# ComplexCruncher
CMPLCRUNCHER=/home/ubuntu/Desktop/complexCruncher/complexCruncher/cmplxcruncher.py

# TODO: 
# Create a new custom project folder inside the "data" folder
# Copy the input xlsx file into your project folder
MYPROJECTFOLDER=/home/ubuntu/Desktop/data/miproyecto
cd $MYPROJECTFOLDER

# Run the ComplexCruncher 
python3 $CMPLCRUNCHER -a pdf -p . > output.log

