function test_it(Q, z0, z1, mu, t, A, b, X0)

constraints_A15_1 = Q; #>= 0

if all(eig(constraints_A15_1) >= 0)
else
  disp("Q NOT pos semidef")
  disp(eig(constraints_A15_1))
end
# Q is Semidefinite

# (A.10)
for i = 1:size(X0,2)
  constraints_A15_2 = [Q,       X0(:,i);
                       X0(:,i)', 1      ]; #>= 0    
  if all(eig(constraints_A15_2) >= 0)
  else
    disp("X0 NOT pos semidef, (A.10)")
    disp(eig(constraints_A15_2))
  end
end

# (A.11)
constraints_A15_3 = [Q,   z0;
                     z0',  1 ]; #>= 0
if all(eig(constraints_A15_3) >= 0)
else
  disp("Q,z0 NOT pos semidef, (A.11)")
  disp(eig(constraints_A15_3))
end

# (A.12)
constraints_A15_4 = [Q,   z0;
                     z0',  mu**2 ]; #>= 0
if all(eig(constraints_A15_4) >= 0)
else
  disp("Q,z1 NOT pos semidef, (A.12)")
  disp(eig(constraints_A15_4))
end


# (A.13)
constraints_A15_5 = Q*A + A'*Q - b*z0' - z0*b'; #<=0
if all(eig(constraints_A15_5) <= 0)
else
  disp("Q,A,z0 NOT neg semidef, (A.13)")
  disp(eig(constraints_A15_5))
end

# (A.14)
constraints_A15_6 = Q*A+A'*Q-b*z0'-z0*b'+2*t*Q; #<=0
if all(eig(constraints_A15_6) <= 0)
else
  disp("Q,A,z1 NOT neg semidef, (A.14)")
  disp(eig(constraints_A15_6))
end
