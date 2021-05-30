format short;
clear all;
clc;
%||| convert to maximization |||
%||| convert all constraints to <= |||

variables = {'x1','x2','x3','s1','s2','sol'}; %includes last as sol
cost = [-2 0 -1 0 0 0]; %cost as in Z include last as sol
info = [-1 -1 1;-1 2 -4]; %coefficients of constraints without slack
b = [-5;-8];  %RHS of constraints
s = eye(size(info,1)); %identity matrix corresponding to slack variables
A = [info s b]; %combine to get the table matrix
BV = [4 5]; %set the basic variables according to identity created in A
fprintf('initial basic variables are : \n');
disp(variables(BV));
%calculate Zj - Cj
ZjCj = cost(BV)*A - cost;
%printing intial table
ZCj = [ ZjCj ; A ];
initialTable = array2table(ZCj);
initialTable.Properties.VariableNames(1:size(A,2)) = variables;
disp(initialTable);

RUN = true;
while RUN
%get last column of  A
sol = A(:,end);
%if any value is negative then not feasible
if any(sol<0)
    fprintf('current BFS is not feasible\nproceeding with iteration\n');
    %finding leaving variable
    [leaveVal pvt_row] = min(sol); %most negative value
   disp('pivot row is : ');
   disp(pvt_row);
    %finding entering variable
    %ratio of pvt_row and Zj - Cj
    row = A(pvt_row,1:end-1); %leave out last column
    ZJ = ZjCj(:,1:end-1); %leave out last column
    for i=1:size(row,2) %for every element in row
        if row(i)<0 %if negative
            ratio(i) = abs(ZJ(i)./row(i)); %take mod of ratio
        else  %if positive assign large value
            ratio(i) = inf;
        end
    end
    [minVal pvt_col] = min(ratio);
    pvt_key = A(pvt_row,pvt_col);
    disp('pivot column is : ');
    disp(pvt_col);
  	BV(pvt_row) = pvt_col;
    disp('basic variables are now : ');
    disp(variables(BV));
    %update the matrix
    A(pvt_row,:) = A(pvt_row,:)./pvt_key; %make pvtkey as 1
    for i=1:size(A,1)       %make all other elements in this column as 0
            if i~=pvt_row
                A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
            end
    end
       %calculate new Zj - Cj
       ZjCj = cost(BV)*A - cost;
       %print this iteration in a table
       ZCj = [ZjCj;A];
       simpTable = array2table(ZCj);
       simpTable.Properties.VariableNames(1:size(ZCj,2)) = variables;
       disp(simpTable);
else
    RUN = false;
    fprintf('This Solution is feasible and optimal\n');
    final_BFS = zeros(1,size(A,2));
    final_BFS(BV) = A(:,end);
    final_BFS(end) = sum(final_BFS.*cost);
    optimalBFS = array2table(final_BFS);
    optimalBFS.Properties.VariableNames(1:size(optimalBFS,2))=variables
end
end