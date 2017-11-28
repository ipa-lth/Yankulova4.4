# A: current system matrix in canonical form ("Steuerung normal form")
A1 = [ -4,      -2;
        1,      -1];
          
 # b: input vector, what does it means?
b1  =   [1; 0];

c1 = [4; 0];     

% [num, den]   = ss2tf(A1, b1, c1, 0);
% [A,b,c,d]    = tf2ss(num,den);
[A,b,c,d,Ti,Qi] = get_Steuerungsnormalform(A1, b1, c1, 0);

%result
%n =
% 2
%Q =
%   1  -4
%   0   1
%Q_inv =
%   1   4
%   0   1
%t1 =
%
%   0   1
%
%T =
%
%   0   1
%
%A0 =
%
%   0   1
%  -6  -5
%
%b0 =
%
%   0
%   1
%
%c0 =
%
%   4
%   4