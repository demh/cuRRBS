#!/bin/bash

###############################################################################
#### Wrapper script for cuRRBS ####
#### Include description ###

#INPUT
# -s: absolute path to the main cuRRBS software folder (e.g. ~/Desktop/cuRRBS)
# -o: absolute path to main output folder (e.g. ~/Desktop/my_cuRRBS_run) 
# 


VERSION=0.1


############################### FUNCTIONS #####################################

###### AUXILIARY FUNCTIONS ######

# print_help()
# Prints a help guide

print_help() {

	echo
	echo
	echo "This is the help page"
	echo
	echo

}



###### PIPELINE STEPS ######



############################## INITIAL MESSAGING ##############################


############################### DEPENDENCIES CHECK ############################

echo
echo
echo "Checking the dependencies ..."
echo

###### PYTHON LIBRARIES ######

PYTHON_D=( os sys pyfaidx Bio csv numpy collections pandas )

for lib in "${PYTHON_D[@]}"
do
	python -c "import $lib" 2>/dev/null
	
	if [ $? -eq 1 ]
	then
	echo "ERROR: The Python library $lib is not installed. Please install it so the software can run."
	echo
	exit 1
	fi 	

done 


###### R PACKAGES ######



echo "The dependencies are correctly installed."
echo

############################### HANDLING ARGUMENTS ############################


############################### RUNNING PIPELINE ##############################



############################ END OF THE SCRIPT ################################

  
