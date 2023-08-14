# InttEventDisplay

How to use INTT Event Display

1. Compere an event base root file which 8felix data conbined
2. Convert to dst file 
    1. Command 
        1. root Fun4All_Anatutorial_Intt.C 
    2. Programs to use
        1. /sphenix/u/mfujiwara/Workspace/tutorials/inttgitclone/AnaTutorial/src/InttAna.h 
        2. /sphenix/u/mfujiwara/Workspace/tutorials/inttgitclone/AnaTutorial/src/InttAna.cc
        3. /sphenix/u/mfujiwara/Workspace/tutorials/inttgitclone/AnaTutorial/macro/Fun4All_Anatutorial_Intt.C
            * This code sets input file by L503, InttRawData *inttraw = new InttRawData("Inttraw", inputfilename);
            * This code sets outputfile name by L59, const string &outputFile = outputfilename,
              
3. Use Event display
    1. Command 
        1. root Loadfile.C 
	2. inttEventDislay->Display_3D()or Display_rphi()or Display_rhoz()
	   Display_3D() shows 3D view and you can turn event display.
	   Display_rphi() shows figure like x-y plane.
	   Display_rhoz() shows figure like y-z plane.
    2. Programs to use
        1.InttEventDisplay/intteventdisplay
            * InttEventDisplay/intteventdisplay/macro/Loadfile.C sets inputfile by L51.

Please refer to this web page to build the package of fun4all
https://github.com/sPHENIX-Collaboration/tutorials/tree/master/AnaTutorial

rootDstDecoder in this repository dosen't work in preparation.
