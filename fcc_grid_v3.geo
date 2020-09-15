SetFactory("OpenCASCADE");
Lx = 18;
Ly = 16;
Mesh.CharacteristicLengthMin = 0.2;
Mesh.CharacteristicLengthMax = 0.2;


Point(1) = {-Lx*0.5, -Ly*0.5, 0};
Point(2) = {Lx*0.5, -Ly*0.5, 0};
Point(3) = {Lx*0.5, Ly*0.5, 0};
Point(4) = {-Lx*0.5, Ly*0.5, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Periodic Line{1} = {-3};
Periodic Line{2} = {-4};

r = 1;
rad=0.85*r;
Disk(2) = {0, 0, 0, rad};
Disk(3) = {2*r, 0, 0, rad};
Disk(4) = {4*r, 0, 0, rad};
Disk(5) = {6*r, 0, 0, rad};
Disk(6) = {8*r, 0, 0, rad};
Disk(7) = {-2*r, 0, 0, rad};
Disk(8) = {-4*r, 0, 0, rad};
Disk(9) = {-6*r, 0, 0, rad};
Disk(10) = {-8*r, 0, 0, rad};

counter=11;

i=1;
For count In {1:2}
    For j In {1:9:2}
        Disk(counter) = {-j*r, -1.732*i*r, 0, rad};
        counter=counter+1;
        Disk(counter) = {j*r, -1.732*i*r, 0, rad};
        counter=counter+1;   
        Disk(counter) = {-j*r, 1.732*i*r, 0, rad};
        counter=counter+1;
        Disk(counter) = {j*r, 1.732*i*r, 0, rad};
        counter=counter+1; 
        Printf("count = %f",count) >>"text.txt";
    EndFor
    i=i+2;
EndFor

i=2;
For count In {1:2}
    Disk(counter) = {0, -1.732*i*r, 0, rad};
    counter=counter+1;
    Disk(counter) = {0, 1.732*i*r, 0, rad};
    counter=counter+1;
    For j In {2:8:2}
        Disk(counter) = {-j*r, -1.732*i*r, 0, rad};
        counter=counter+1;
        Disk(counter) = {j*r, -1.732*i*r, 0, rad};
        counter=counter+1;   
        Disk(counter) = {-j*r, 1.732*i*r, 0, rad};
        counter=counter+1;
        Disk(counter) = {j*r, 1.732*i*r, 0, rad};
        counter=counter+1; 
        Printf("count = %f",count) >>"text.txt";
    EndFor
    i=i+2;
EndFor



f1() = BooleanDifference{ Surface{1}; }{ Surface{2:86}; };
f2() = BooleanIntersection{ Surface{1}; Delete; }{ Surface{2:86}; Delete; };

Physical Surface(0)={2,3,7,10,12,13,15,16,19,20,21,25,26,32,33,35,36,39,40,41,45,46,51,53,54,60,62,63,65,66,69,70,73,76,77,81,82,83,85};
Physical Surface(1) = {87};
Physical Surface(2) = {4,5,6,8,9,11,14,17,18,22,23,24,27,28,29,30,31,34,37,38,42,43,44,47,48,49,50,52,55,56,57,58,59,61,64,67,68,71,72,74,75,78,79,80,84,86};

Mesh 3;
Coherence Mesh;
Mesh 2;
