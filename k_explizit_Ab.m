function [k0, k1] = k_explizit_Ab(roots_p1, roots_pmin, pmin, A, b)

[A_R, _, _, _ , _, _] = get_Steuerungsnormalform(A, b, b, 0);
a = -A_R(size(A_R, 1), :)';
n = size(roots_p1, 2);

a_tilde_p1 = fliplr(poly(roots_p1)(2:n+1))'; # get the characteristical polynomial backwards
a_tilde_pmin = fliplr(poly(roots_pmin)(2:n+1))'; # get the characteristical polynomial backwards

k1 = (a_tilde_p1 - a_tilde_pmin) / (1.0-1.0/pmin);
k0 = a_tilde_p1 - a - k1;

return