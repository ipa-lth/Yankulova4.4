function [R_sol, sol, res] = solveLMIs(inStruct)

%%
% input parameters
n = inStruct.n;
p_min = inStruct.p_min;
bet = inStruct.bet;
M_l = inStruct.M_l;
M_u = inStruct.M_u;
N_l = inStruct.N_l;
N_u = inStruct.N_u;
% ScaleFactor = inStruct.ScaleFactor;

% plant and contoller
for i1 = 0:(-N_l+N_u)
    k{i1+1} = inStruct.k{i1+1}; %#ok<AGROW>
    A{i1+1} = -inStruct.B*k{i1+1}'; %#ok<AGROW>
    ks{i1+1} = inStruct.ks{i1+1}; %#ok<AGROW>
    As{i1+1} = -inStruct.B*ks{i1+1}'; %#ok<AGROW>
end
A{-N_l+1} = A{-N_l+1} + inStruct.A;  % term with p^0!
As{-N_l+1} = As{-N_l+1} + inStruct.A;  % term with p^0!
if N_u<0
    error('N_u<0 not yet implemented!')
end
x0 = 1.001*inStruct.x0;


%% define and solve LMI-Problem

% Q_i, the actual lmi-variables
for i1 = 0:(M_u-M_l)
    R{i1+1} = sdpvar(n,n,'symmetric'); %#ok<AGROW>
end

% % condition a) R_{M_l} > 0
% LMIs = set(R{1} > 0);

% transform LMI variables for checking condition b)
% Q(p) positive definite for all p \in [p_min, 1]
R_b = transQ(R, (1-p_min)/2, (1+p_min)/2);
R_b_sum = -calcQsum(R_b); % negative since constraint is positive definite
lmi_out = calcLMIcond(R_b_sum,n);
LMIs = lmi_out;

% transform LMI variables for checking condition d)
% nesting condition II
for i1 = 0:(-M_l+M_u)
    R_d{i1+1} = (i1+M_l) * R{i1+1}; %#ok<AGROW>
end
R_d = transQ(R_d, (1-p_min)/2, (1+p_min)/2);
R_d_sum = calcQsum(R_d);
lmi_out = calcLMIcond(R_d_sum,n);
LMIs = [LMIs, lmi_out];

% transform LMI variables for checking condition e1)
% stability condition
for j1 = 0:(-M_l+M_u-N_l+N_u)
    R_e{j1+1} = 0; %#ok<AGROW>
    for i1 = 0:(-M_l+M_u)
        i2 = j1-i1;
        if (i2>=0)&&(i2<=(-N_l+N_u))
            R_e{j1+1} = R_e{j1+1} + R{i1+1}*A{i2+1} + A{i2+1}'*R{i1+1}; %#ok<AGROW>
        end
    end
end
R_e = transQ(R_e, (1-p_min)/2, (1+p_min)/2);
R_e_sum = calcQsum(R_e);
lmi_out = calcLMIcond(R_e_sum,n);
LMIs = [LMIs, lmi_out];

% transform LMI variables for checking condition e2)
% stability condition (saturated system)
for j1 = 0:(-M_l+M_u-N_l+N_u)
    R_e2{j1+1} = 0; %#ok<AGROW>
    for i1 = 0:(-M_l+M_u)
        i2 = j1-i1;
        if (i2>=0)&&(i2<=(-N_l+N_u))
            R_e2{j1+1} = R_e2{j1+1} + R{i1+1}*As{i2+1} + As{i2+1}'*R{i1+1}; %#ok<AGROW>
        end
    end
end
R_e2 = transQ(R_e2, (1-p_min)/2, (1+p_min)/2);
R_e2_sum = calcQsum(R_e2);
lmi_out = calcLMIcond(R_e2_sum,n);
LMIs = [LMIs, lmi_out];

% condition f1)
% actuator constraint
L = min(N_l, M_l);
for j1 = 0:(max(N_u, M_u) - min(0, L))
    R_f{j1+1} = 0; %#ok<AGROW>
    if (j1<=M_u-L) && (j1>=M_l-L)
        R_f{j1+1} = R_f{j1+1} + [[0, zeros(1,n)]; [zeros(n,1), [R{j1-(M_l-L)+1}]]];
    end
    if (j1<=N_u-L) && (j1>=N_l-L)
        R_f{j1+1} = R_f{j1+1} + [[0, k{j1-(N_l-L)+1}']; [k{j1-(N_l-L)+1}, zeros(n)]];
    end
end
R_f{-min(0, L)+1} = R_f{-min(0, L)+1} + [1, zeros(1,n); zeros(n,1), zeros(n)];
R_f = transQ(R_f, (1-p_min)/2, (1+p_min)/2);
R_f_sum = -calcQsum(R_f); % negative since constraint is positive definite
lmi_out = calcLMIcond(R_f_sum,n+1);
LMIs = [LMIs, lmi_out]; 

% condition f2)
% actuator constraint
R_1 = R{1};
for i1 = 1:(-M_l+M_u)
    R_1 = R_1 + R{i1+1};
end
ks_1 = ks{1};
for i1 = 1:(-N_l+N_u)
    ks_1 = ks_1 + ks{i1};
end
LMIs = [LMIs, [[[bet^2, ks_1']; [ks_1, R_1]] > 0]];

% condition g)
% initial states constraint
for i1 = 1:size(x0,2)
    LMIs = [LMIs, [1 - x0(:,i1)'*R_1*x0(:,i1) > 0]];
end

% objective function
R_min = 0;
for i1 = M_l:M_u
    R_min = R_min + p_min^i1*R{i1+1-M_l};
end
    
% for i1 = 1:size(x0,2)
%     LMIs = LMIs + set(1 - ScaleFactor^2*x0(:,i1)'*R_min*x0(:,i1) > 0);
% end
% h = -geomean(R_min);
h = trace(R_min);
% h = diag(R_min)'*diag(R_min);

opts = sdpsettings('solver','sdpt3','sdpt3.gaptol',1e-5,'verbose',0, 'savesolveroutput', 1,'cachesolvers',1);
% opts = sdpsettings('solver','sedumi','sedumi.numtol',1e-7, 'savesolveroutput', 1);
% opts = sdpsettings('solver','csdp', 'savesolveroutput', 1);
% opts = sdpsettings('solver','dsdp','dsdp.gaptol',1e-7);
sol = solvesdp(LMIs, h, opts);
for i1 = 0:(-M_l+M_u)
    R_sol{i1+1} = double(R{i1+1}); %#ok<AGROW>
end
[pres,dres] = checkset(LMIs);
res = min(min(pres),min(dres))/max(max(pres),max(dres));



end % main 



%% ----------------------
%% ----------------------

function Q_out = transQ(Q_in,a,b)
% calculate the interval transformation (p(q) = a*q+b)

if ~iscell(Q_in)
    error('cell array input expected')
end

m = length(Q_in);
n = size(Q_in{1},1);

for i1 = 0:m-1
    Q_out{i1+1} = Q_in{i1+1}; %#ok<AGROW>
    for i2 = i1+1:m-1
        Q_out{i1+1} = Q_out{i1+1} + b^(i2-i1)*nchoosek(i2,i1)*Q_in{i2+1}; %#ok<AGROW>
    end
    Q_out{i1+1} = a^i1*Q_out{i1+1}; %#ok<AGROW>
end

end % transQ



%% ----------------------

function Q_sum = calcQsum(Q_in)
% calculate the 'blown up' matrix

if ~iscell(Q_in)
    error('cell array input expected')
end

m = length(Q_in);
n = size(Q_in{1},1);
if m == 1
    Q_sum = [2*Q_in{1}, zeros(n); zeros(n), zeros(n)];
elseif m == 2
    Q_sum = [[[2*Q_in{1}], [Q_in{2}]]; [[Q_in{2}], zeros(n)]]
else
    Q_sum = [[[2*Q_in{1}], [Q_in{2}]]; [[Q_in{2}], [2*Q_in{3}]]];
end
for i1 = 4:2:m
    if i1 ~= m
        Q_sum = [[Q_sum, [zeros((i1/2-1)*n,n); Q_in{i1}]]; ...
            [[zeros(n,(i1/2-1)*n), Q_in{i1}], [2*Q_in{i1+1}]]]; %#ok<AGROW>
    else % even number of elementary matrices --> zeros in southeasternmost position
        Q_sum = [[Q_sum, [zeros((i1/2-1)*n,n); Q_in{i1}]]; ...
            [[zeros(n,(i1/2-1)*n), Q_in{i1}], zeros(n)]]; %#ok<AGROW>
    end
end
Q_sum = 1/2*Q_sum;

end % calcQsum


%% ----------------------

function lmi_out = calcLMIcond(P_sum,n)
% calculate lmi condition

k = size(P_sum,1)/n;

J = [zeros(n*(k-1),n), eye(n*(k-1))];
C = [eye(n*(k-1)), zeros(n*(k-1),n)];
D = sdpvar(n*(k-1),n*(k-1),'symmetric');
G = sdpvar(n*(k-1),n*(k-1),'skew');
lmi_out = [D > 0];
lmi_out = [lmi_out, P_sum < [C;J]'*[-D, G; G' D]*[C;J]];

end 
