function Keff = solve_cell_problem(j,user_data,KA,KB)
% ------------------------------------------------------------------------------
% solve_cell_problem.m
%
% Written by Dr Elliot Carr (2013-2014)
% Ecole Centrale Paris and Queensland University of Technology
% 
% This code is part of TwoScalRich.
%
% This MATLAB function solves the periodic cell problem:
%
%            \nabla_y \cdot [ K(h,y)\nabla_y (u_{j} + y_{j}) ] = 0
%
% where K(h,y) = KA in soil A and K(h,y) = KB in soil B, sub to the constraint 
% that the cell-averaged value of u is equal to 0 over the cell.
% ------------------------------------------------------------------------------

m_no_variables   = user_data.micro_mesh.no_variables;
user_data.KA     = KA;
user_data.KB     = KB;
user_data.column = j;

% Problem is linear and has the form F = Au + b = 0, where A = [A1; A2] and 
% b = [F1; F2].
u  = zeros(m_no_variables,1); epsilon = 1.0;
F1 = feval('Ffunc1',u,user_data);

% Build A1
A1 = zeros(m_no_variables,m_no_variables);
e  = zeros(m_no_variables,1);
for i = 1:m_no_variables
    e(i)  = 1;
    A1(:,i) = (feval('Ffunc1',u+epsilon*e,user_data) - F1) / epsilon;
    e(i)  = 0;
end

F2 = feval('Ffunc2',u,user_data);

% Build A2
A2 = zeros(1,m_no_variables);
e = zeros(m_no_variables,1);
for i = 1:m_no_variables
    e(i)  = 1;
    A2(i) = (feval('Ffunc2',u+epsilon*e,user_data) - F2) / epsilon;
    e(i)  = 0;
end
A = [A1; A2];
b = [F1; F2];

% Solution
u = A\(-b);

% Compute column j of tensor Keff(h)
Keff = compute_Keff(u,user_data);

end