

PABLO - POTENTIAL FLOW AROUND AIRFOILS WITH                    /pablo/readme.txt
   BOUNDARY LAYER COUPLED ONE-WAY
       
The files for this program are in the directory pablo 
  readme.txt      this file - general description
  pablo.txt       instructions for usage
  pablo.tar       tar file as downloaded from KTH web site


The following files make up the source and program code:
cfturb.m      close_te.m    dist.m        doublet.m     ellipse.m   
fh.m          findy.m       fl.m          h1ofh.m       hofh1.m     
info.m        library.m     naca4.m       pablo.m       rsolver.m   
runge.m       solvebl.m     source.m      splf.m        stag.m  
sy.m          vortex.m    

The following example cases are included:
sc21006.dat 	clarky.dat  	dae11.dat   	dae21.dat   	dae31.dat   
dsma523.dat 	epp662.dat  	epp748.dat  	foil31.dat  	fx63137.dat 
gaw1.dat    	hsn0213.dat 	k720616.dat 	k790312.dat 	k820609.dat 
korn.dat    	liss7769.dat	ls10013.dat 	ls10417m.dat	ms10313.dat 
ms10317.dat 	n001035.dat 	n63215.dat  	n643418.dat 	n64a010.dat 
n64a410.dat 	n651012.dat 	n651213.dat 	n652215.dat 	n65410.dat  
n658299m.dat	n658299r.dat	n65a012.dat 	n663018.dat 	n747a415.dat
naca0012.dat	naca4412.dat	newabb.dat  	nl10414f.dat	nl10416.dat 
nl11215f.dat	nl20415.dat 	nlr1.dat    	oneram6.dat 	rae100.dat  
rae101.dat  	rae102.dat  	rae103.dat  	rae104.dat  	rae2822.dat 
sc20010.dat 	sc20012.dat 	sc20406.dat 	sc20410.dat 	sc20412.dat 
sc20414.dat 	sc20518.dat 	sc20606.dat 	sc20610.dat 	sc20612.dat 
sc20614.dat 	sc20706.dat 	sc20712.dat 	sc20714.dat 	cast7.dat   
sc21010.dat 	super11.dat 	vezbl32.dat 	vezcan.dat  	vezwltr.dat 
wilbyb.dat  	wilbyc.dat  	wilbyr.dat  
These example files are compressed in the archive file cases.zip
 

Pablo is one of the first contributions to this CD-ROM that did not originate
at one of the U.S. national laboratories. This is a contribution from KTH -
The Royal Institute of Technology, Department of Aeronautics, Stockholm, Sweden.
Pablo was programmed by Christian Wauquiez in 1999. His KTH mentor was
Prof. Arthur Rizzi.

Pablo is a pedagogical low-speed airfoil analysis program written in MATLAB. 
It uses a one way coupled inviscid + boundary layer model. The inviscid 
flow is solved using a Panel Method [1]. Three different kinds of 
singularity distributions can be used. The boundary layer equations use 
the inviscid flow velocity provided by the panel method, but the effect of 
the boundary layer on the inviscid flow is not taken into account, as
in Panda [2]. The boundary layer model is described in [3]. It uses Thwaites'
equations for the laminar part of the flow, and Head's equations for the
turbulent part. Michel's criterion is used to locate transition.
The drag coefficient is computed using the Squire-Young formula.

This program is included on this CD-ROM as a convenience. The definitive
source for the program is the web site
  https://www.nada.kth.se/~chris/pablo/pablo.html
and the program may be downloaded from there.

Pablo is the first program in this collection that was created with
MatLab. Programs written with MatLab are usually quite understandable
and readable because many of the numerical details are hidden from
view. The only drawback is that you need to have MatLab installed on
your computer in order to run the program. You do not get an *.exe
file that can be run anywhere. The good news is that you only need
the student edition of MatLab and not the expensive full edition.

REFERENCES
[1] Katz and Plotkin : Low Speed Aerodynamics, From Wing
Theory To Panel Methods. McGraw-Hill Inc., 1991 
[2] Ilan Kroo : PANDA - A Program for Analysis and Design of
Airfoils. Desktop Aeronautics, 1988.
[3] Jack Moran : An introduction to Theoretical and Computational
Aerodynamics. John Wiley and sons, 1984.
[4] I.H. Abbott and A.E. Von Doenhoff, Theory of wing sections,
Dover Publications Inc., New York, 1959 



KTH- The Royal Institute of Technology
Department of Aeronautics
Stockholm, Sweden

Programmed by Christian Wauquiez, 1999
KTH Contact : Prof. Arthur Rizzi, rizzi@flyg.kth.se
