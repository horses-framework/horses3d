Mesh.MshFileVersion = 4.1;

// ----------- BEGIN INPUT ------------
N = 16;
Nx = 5;
Na = 30;
scale = 1.0;     // geometry scaling
// ----------- END INPUT --------------

l1=1.0;
l2=10.0;
Lx=1.0;
c1 = 1.0;

Point(1) = {0.0, 0.0, 0.0, c1};
Point(2) = {0.0, l1, 0.0, c1};
/* Point(3) = {0.0, l2, 0.0, c1}; */

line[] = Extrude {0.0, l2, 0.0} {
  Point{2}; Layers{N};
};

surface[] = Extrude {Lx, 0.0, 0.0} {
  Line{line[1]}; Layers{Nx}; Recombine;
};
//+
Extrude {{1, 0, 0}, {0, 0, 0}, 2*Pi/3} {
  Surface{5}; Layers{Na}; Recombine;
}
//+
Physical Surface("inflow") = {14};
Physical Surface("outflow") = {22};
Physical Surface("top") = {18};
Physical Surface("bottom") = {26};
Physical Surface("sidea") = {5};
Physical Surface("sideb") = {27};
//+
Physical Volume("fluid") = {1};
//+
Mesh.ElementOrder = 3;
