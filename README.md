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
    1. Set data file.
   	InttEventDisplay/macro/Loadfile.C set inputfile name by L51. 
    3. Enter the Commands in the order.
        1. root Loadfile.C
   		Data loaded from the begining and stop run when it found a event, Ncluster > 2.
	4. inttEventDislay->display()
    		Display x-y and y-z plane and 3D viewer. 3DViewer is displayed another tab.
    	3. se->run(0)
   		Resume data load when find a next event, Ncluster > 2.
