

LINEINT - INTERSECTION OF TWO STRAIGHT LINES        /lineint/readme.txt

The files for this program are in the directory lineint 
  readme.txt      general description
  xsect.dfm       the main form
  xsect.pas       the code for unit main 
  int2d.dpr       the Delphi project file
  int2d.exe       the executable file for Windows
  xlines.c        original article from Graphics Gems (not used)
  xsect.html      web page with JavaScript calculation of intersection 

In order to recompile this program, you will need to have Delphi installed
on your computer. You can use the executable file int2d.exe without Delphi.

Int2D is a simple geometry tool to calculate the intersection of two
lines. Each line is defined by two points. This is hardly a difficult
math problem, but you can eat up a lot of time and make errors working
out a set of intersections with pencil and calculator. I find this
program quite useful in doing layouts of various concepts and components.

The code is based on an article in the book Graphics Gems, describing the
best way to calculate the intersection - with the fewest number of operations.
This turns out to be 10 multiplications, 9 subtractions and 2 divisions.
Can you do better?? You might find the code in xsect.pas interesting.
You can pull out the procedure called LinesIntersect and use it if
you are writing a program and need to compute an intersection.

The usage is supposed to be self-evident. Fill in the blanks and
press "Compute". When you are done, press Quit. You get a simple
error message if you input parallel lines or try to define a line with
two identical points. 

This program is written in Pascal (Delphi) and only has an executable 
for Windows. 

I wrote this program as a new project and is not based on a program
from a national laboratory. The whole program, source code and all is
declared to be open source - public domain.

For October 2020, I have added the algorithm as a self-contained website
that includes the calculation as a web script (JavaScript).
This page is also reachable via the PDAS website as 
      https://www.pdas.com/xsect.html
You can reach this page from a smartphone or tablet or a computer of any
operating system, thereby making the program universally available.      
You may not know this, but if you have xsect.html stored on your machine,
then you launch the web page by double clicking on the file name. 
This works even if you are not currently connected to the internet.
