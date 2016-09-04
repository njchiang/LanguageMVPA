#!/bin/bash
:<<doc
Written 10/7/2014 by Jeff Chiang

this function takes the full path of a matlab script and runs it on the grid. eventually it will also take parameters
note: script has to have all subroutines and directories defined (setpath and addtopath), and executable permissions


matlab script will have to have the following lines at the top, this is to avoid any path issues
------------------
toolboxRoot=[PATH_WHERE_SCRIPTS_ARE]; 
addpath(genpath(toolboxRoot))
cd [DIRECTORY OF CHOICE]

##INITIALIZE YOUR PARAMETERS
##CALL YOUR FUNCTION WITH PARAMETERS... they have to be hard-coded for now but it's not too much work to modify a .m script =) 
------------------

to take full advantage of this, use PARFOR. But I won't write a guide on that unless I have time.

doc

scriptName=$1
currDir=/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/grid
cd ${currDir}

echo starting script ${scriptName}

echo "Writing temporary script..."
cat <<EOF > temp.sh  
#!/bin/bash 
matlab -nodesktop -nosplash -nodisplay -r "run ${scriptName} ; quit"
EOF
chmod 775 temp.sh
sge qsub ${currDir}/temp.sh
echo done submitting ${scriptName} 
echo "remove temp.sh?"
rm -i temp.sh
 
