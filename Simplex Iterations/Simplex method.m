format short;
clear all;
clc;
%||| simplex method |||
%||| all constraints should be <= |||
%||| problem should be maximization accordingly adjust C |||

C = [-1,3,-2]; % vector for Z
noVariables = 3;
info = [3,-1,2;-2,4,0;-4, 3,8]; % coefficients in constraints
b = [7;12;10] ; %RHS of constraints
%identity matrix of the size equal to number of constraints
%size(info,1) gives number of rows in info matrix
s = eye(size(info,1)) ;
% forming the starting simplex table
A = [info s b];
%cost matrix for the table (Cj matrix)
Cost = zeros(1,size(A,2)); %size(A,2) defines the number of columns in A
Cost(1:noVariables) = C; %setting the cost for variables and setting the cost of slack variables as 0
%initial basic variables as slack variables
BV = noVariables+1:1:size(A,2)-1;  %column index after starting variables and leaving last column for b

%calculating Zj - Cj (costs of basic variables being multiplied with each
%column and subtracting Cj
ZjCj = Cost(BV)*A -Cost; %Cost(BV) implies the cost vector for column numbers defined in BV

%printing the complete table
ZCj = [ZjCj;A];  %A matrix below Zj - Cj row
simpTable = array2table(ZCj); %better representation
simpTable.Properties.VariableNames(1 : size(ZCj,2)) = {'x1','x2','x3','s1','s2','s3','sol'};
simpTable

%simplex method starts
RUN = true;
while RUN

%creating further optimal tables
if any(ZjCj<0)
    disp('not optimal table');
    disp('old basic variables : ');
    disp(BV);
    
    %finding pivot column using min ZjCj value
    ZC = ZjCj(1 : end-1);
    [EnterCol,pvt_col] = min(ZC) %EnterCol stores the minimum ZjCj value and pvt_col stores the index of this column
    fprintf('min Zj - Cj value is : %d ',EnterCol);
    fprintf('\nentering variable is : x%d \n',pvt_col);
    
    %finding minimum ratio using b/xi as Q for pivot row
    sol = A(:,end);   %last column of the table as b A(:,end) implies all row and end column
    column  = A(:,pvt_col);  %column for entering variable
    %if all values in column are negative then unbounded solution
    if all(column<0)
        error('LPP is unbounded all entries are zero in column : %d',pvt_col);
    end
    %calculate ratio only for positive values
    for i=1:size(column,1)
        if column(i) > 0
            ratio(i) = sol(i)./column(i);   % ./ is to divide each element with corresponding element
        else
            ratio(i) = inf;  %assign invalid value as inf easier to choose min ratio then
        end
    end
    %find minimum ratio
    [minRatio , pvt_row] = min(ratio);
    fprintf('min ratio is : %d \npvt_row is : %d \n',minRatio,pvt_row);
    fprintf('leaving variable is : %d\n',BV(pvt_row));
    %change basic variables according to leaving and entering and select
    %pvt_element
    BV(pvt_row) = pvt_col;
    disp('new basic variables : ');
    disp(BV);
    pvt_key =  A(pvt_row,pvt_col);
    
    %update the table
    A(pvt_row,:) = A(pvt_row,:)./pvt_key;  %row operation making pivot element as 1
    %row operation making all other elements in pvt_col as 0
    for i=1:size(A,1)
        if i~=pvt_row
            A(i,:) = A(i,:) - A(i,pvt_col).*A(pvt_row,:);
        end
        ZjCj = ZjCj - ZjCj(pvt_col).*A(pvt_row,:);  %same operation for Zj - Cj row
    end
    %printing the complete table
    ZCj = [ZjCj;A];  %A matrix below Zj - Cj row
    simpTable = array2table(ZCj); %better representation
    simpTable.Properties.VariableNames(1 : size(ZCj,2)) = {'x1','x2','x3','s1','s2','s3','sol'};
    simpTable
    %print the current BFS
    BFS = zeros(1,size(A,2));
    BFS(BV) = A(:,end);
    BFS(end)  = sum(BFS.*Cost);
    current_BFS = array2table(BFS);
    current_BFS.Properties.VariableNames(1:size(current_BFS,2)) = {'x1','x2','x3','s1','s2','s3','sol'};
    current_BFS
    %if problem is minimization then print the negative of the final BFS
    %sol value
    
else
    RUN=false;
    disp('optimal solution reached');
end
end
