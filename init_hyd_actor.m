###########################
# Hydraulischer Aktor     #
###########################

A0 = [    0,        1,       0;
        -10,   -1.167,      25;
          0,        0,    -0.8];

b0 = [0; 0; 2.4];
c0 = [1; 0; 0];
d0 = [0];

u_max = 10.5;
n = 3;

X00     = [ -20, -10, -10;
            -20, -10,  10;
            -20,  10, -10;
            -20,  10,  10;
             20, -10, -10;
             20, -10,  10;
             20,  10, -10;
             20,  10,  10];

# A, b,c  should be in canonical form ("Steuerung normal form")
[A, b, c, d, T, Q] = get_Steuerungsnormalform(A0, b0, c0, d0);

a = (-1)*A(size(A)(2),:)'; # a.T is the last line of canonical form of A 

X0 = (T*(X00.'));
