% Example: min |Ax-b|

% setup several equations in one dimension which are inconsistent.


A = [1 -1 3 -2];
b = -[5 7  1 1];

% plot the lines 
x = linspace(-3,3); 
xo = ones(size(x));
f = x'*A - xo'*b;
plot(x, f,x, xo'*0,'k');

%% find the solution in different norms
fprintf(1,'Computing optimal solution in the case of Linfty-norm...');

cvx_begin quiet
    variable x1
    minimize (norm(A'*x1 - b', inf))
cvx_end

fprintf(1,'Done! \n');


%%
y1 = A'*x1-b';
figure(1), plot(x, f,x, xo'*0,'k'); 
hold on, plot(x1,y1,'*k'), title('shows error'), 

%%
fprintf(1,'Computing optimal solution in the case of L2-norm...');

cvx_begin quiet
    variable x1
    minimize (norm(A'*x1 - b', 2))
cvx_end

fprintf(1,'Done! \n');

y2 = A'*x1-b';
figure(1), plot(x, f,x, xo'*0,'k'); 
hold on, plot(x1,y2,'*r'), title('shows error'), hold off
norm(y1)
norm(y2)

