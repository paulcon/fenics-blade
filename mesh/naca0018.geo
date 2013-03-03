Function NACA4

// generate a symmetric NACA four series aerofoil composed of a set of
// splines
//
// the following variables must be set on entry:
// NACA4_th     thickness in percent of chord
// NACA4_ch     aerofoil chord
// NACA4_le_x,y,z   leading edge coordinates
// NACA4_len_te length scale (trailing edge)
// NACA4_len_mc length scale (mid chord)
// NACA4_len_le length scale (leading edge)
// NACA4_nspl   number of splines on section
// NACA4_pps    number of points per spline

// The local scale length will be set using a quadratic which interpolates
// trailing edge, midpoint and leading edge scale lengths as a function of
// distance along the chord.

// On exit, the following variables will contain the details of
// a closed NACA 00TH section with 2*NACA4_nspl splines, each of
// NACA4_pps points:
//
// NACA4_Points[] a list of the 2*NACA4_nspl*NACA4_pps points
//                on the section
// NACA4_Splines[] a list of the 2*NACA4_nspl splines
//
// these two lists are oriented so that they start at the trailing edge
// and move over the upper surface and around the lower surface to return
// to the trailing edge

// constants from NASA TM4741
NACA4_a0 =  0.2969 ;
NACA4_a1 = -0.1260 ;
NACA4_a2 = -0.3516 ;
NACA4_a3 =  0.2843 ;
NACA4_a4 = -0.1015 ;

NACA4_npts = NACA4_nspl*NACA4_pps ;

NACA4_Points[] = {} ;
NACA4_Splines[] = {} ;

For i In {0:NACA4_npts}
    x = 0.5*(1+Cos(i/NACA4_npts*Pi)) ;
    y = NACA4_a0*Sqrt(x) + 
        x*(NACA4_a1 + x*(NACA4_a2 + x*(NACA4_a3 + x*NACA4_a4))) ;
    p = newp ;
    L1 = 1-x ; L2 = x ;
    NACA4_len = NACA4_len_le*L1*(2*L1-1) + NACA4_len_mp*4*L1*L2 + NACA4_len_te*L2*(2*L2-1) ;
    Point(p) = {NACA4_ch*x+NACA4_le_x,  
                y*NACA4_ch*NACA4_th/20+NACA4_le_y, NACA4_le_z, NACA4_len} ;
    NACA4_Points[i] = p ;
EndFor

For i In {NACA4_npts+1:2*NACA4_npts-1}
    x = 0.5*(1+Cos(i/NACA4_npts*Pi)) ;
    y = NACA4_a0*Sqrt(x) + 
        x*(NACA4_a1 + x*(NACA4_a2 + x*(NACA4_a3 + x*NACA4_a4))) ;
    p = newp ; 
    NACA4_len = NACA4_len_le*L1*(2*L1-1) + NACA4_len_mp*4*L1*L2 + NACA4_len_te*L2*(2*L2-1) ;
    Point(p) = {NACA4_ch*x+NACA4_le_x,  
                -y*NACA4_ch*NACA4_th/20+NACA4_le_y, NACA4_le_z, NACA4_len} ;
    NACA4_Points[i] = p ;
EndFor

For i In {0:2*NACA4_nspl+NACA4_pps-1}
    c = newc ;
    Spline(c) = NACA4_Points[{i*(NACA4_pps-1):(i+1)*(NACA4_pps-1)}] ;
    NACA4_Splines[i] = c ;
EndFor
c = newc ; i = 2*NACA4_nspl+NACA4_pps ;

Spline(c) = {NACA4_Points[{i*(NACA4_pps-1):(2*NACA4_npts-1)}], NACA4_Points[0]} ;

NACA4_Splines[i] = c ;

Return


NACA4_nspl = 8 ;
NACA4_pps = 4 ;

NACA4_len_le = 0.1 ;
NACA4_len_mp = 0.2 ;
NACA4_len_te = 0.1 ;
NACA4_th = 18 ;
NACA4_ch = 1.0 ;
NACA4_le_x = 0.0 ;
NACA4_le_y = 0.0 ;
NACA4_le_z = 0.0 ;

Call NACA4 ;

Point(65) = {0.2, 0.0, 0, 1.0};
Point(66) = {0.45, 0.0, 0, 1.0};
Point(67) = {0.7, 0.0, 0, 1.0};

Point(68) = {0.15, 0.0, 0, 1.0};
Point(69) = {0.4, 0.0, 0, 1.0};
Point(70) = {0.67, 0.0, 0, 1.0};

Point(71) = {0.25, 0.0, 0, 1.0};
Point(72) = {0.5, 0.0, 0, 1.0};
Point(73) = {0.73, 0.0, 0, 1.0};

Point(74) = {0.2, 0.05, 0, 1.0};
Point(75) = {0.45, 0.05, 0, 1.0};
Point(76) = {0.7, 0.03, 0, 1.0};

Point(77) = {0.2, -0.05, 0, 1.0};
Point(78) = {0.45, -0.05, 0, 1.0};
Point(79) = {0.7, -0.03, 0, 1.0};

Circle(24) = {74, 65, 68};
Circle(25) = {68, 65, 77};
Circle(26) = {77, 65, 71};
Circle(27) = {71, 65, 74};

Circle(28) = {75, 66, 72};
Circle(30) = {72, 66, 78};
Circle(31) = {78, 66, 69};
Circle(32) = {69, 66, 75};

Circle(33) = {70, 67, 76};
Circle(34) = {76, 67, 73};
Circle(35) = {73, 67, 79};
Circle(36) = {79, 67, 70};
Line Loop(37) = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 1, 2, 3, 4, 5};
Line Loop(38) = {27, 24, 25, 26};
Line Loop(39) = {32, 28, 30, 31};
Line Loop(40) = {33, 34, 35, 36};
Plane Surface(41) = {37, 38, 39, 40};

// Extrude {0, 0, 2} {
//   Surface{41};
//   Layers{40};
//   Recombine;
// }

Extrude {0, 0, 2} {
  Surface{41};
  Layers{40};
}

Physical Surface(209) = {208};
Physical Surface(210) = {41};
Physical Surface(211) = {115, 123, 119, 131, 127, 135, 139, 147, 143, 159, 151, 155, 79, 83, 107, 87, 91, 99, 103, 111, 95};
Physical Surface(212) = {171, 167, 175, 163};
Physical Surface(213) = {179, 191, 183, 187};
Physical Surface(214) = {207, 203, 195, 199};
Physical Volume(215) = {1};
