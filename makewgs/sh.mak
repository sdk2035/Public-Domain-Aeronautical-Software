--- This example makes a body with a nose 20 units long and an afterbody
--- that is 10 units in length. Each is shaped as a Sears-Haack body.
--- There will be 21 stations defining the body and 9 meridian lines.
--- on the right half of the body. If used as input to PanAir, there
--- would be 20x8 or 160 panels.

NAME BODY-1
LNOSE 20
LTAIL 10
LBODY 30   !  since LBODY-LNOSE-LTAIL=0 there will be no cylinder portion
RADIUS 3   ! since RNOSE and RBASE are not entered, they are zero
NFS 21
NTHETA 9
BODY
