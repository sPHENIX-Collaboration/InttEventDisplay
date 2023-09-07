# InttEventDisplay

How to use INTT Event Display

1. Compere an event base root file which 8felix data conbined  
2. Convert to dst file  
3. Make build directory and run autogen.sh.  
	1. mkdir build  
 	2. cd build 
 	3. ../autogen --prefix=$MYINSTALL  
	MYINSTALL is your install directory. You can set this enviroment variable by export function in bash shell.  
  	4. make && make install  
  
       	Please refer to this web page to build the package of fun4all  
	https://github.com/sPHENIX-Collaboration/tutorials/tree/master/AnaTutorial

4. Use Event display
 Enter the Commands in the order.
   1. bash InttEventDisplay.sh filepath NClusters true or false 
	DST file is read from the beginning until an event with more clusters than NCluster is found, and then an event display window opens. 
	Display x-y and y-z plane and 3D viewer. 3DViewer is displayed another tab. 
	If the 3rd argument is true, the image is saved; if false, it is not saved. 
   2. se->run(0) 
      	Resume data load until find a next event, which has clusters bigger than NClusters.
