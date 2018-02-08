i = 1;

Q = [1, 1; 1, 1];

a = [1; 1];

n = 2;

I = eye(n);

N = [-2, 0; 0, -1];

M = [0, 0; 0, 1];

z = [1; 1];

tmp = func_constraint464(i, Q, a, M, N, n, z);
disp(tmp);
if (tmp == 9)
  disp("RIGHT calculated!")
else
  disp("WRONG calculated!")
endif