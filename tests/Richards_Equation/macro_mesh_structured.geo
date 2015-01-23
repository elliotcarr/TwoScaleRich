a = 1.0;
b = 1.0;
elements_x = 10;
elements_y = 10;
left = 0;
right = 0.4;
r = a/elements_x;
Point(1) = {0,0,0,r};
Point(2) = {a,0,0,r};
Line(1) = {1,2};
Extrude {0,b,0} {Line{1}; Layers{elements_y}; Recombine; }
