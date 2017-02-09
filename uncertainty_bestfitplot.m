%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%        Bestfit lines
%        Uncertainty
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cleaning up and setup formating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
format short

% Create a figure to plot what this data look like
figure (1)
hold on;
set (gcf,'defaultaxesfontsize',14)
title ('Best fit comparison','FontName','Arial','FontSize',20)
xlabel('time','FontName','Times','FontSize',18,'FontWeight','bold')
ylabel('Y-axis','FontName','Times','FontSize',18,'FontWeight','bold')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% values of x, and y we need to least mean square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = [4,3,1];
t = [-2,-1,0];
[~,c]=size(t);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formulas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%formulas
bestfit=@(j) (AB_val(3)+(AB_val(2)*j));
A_val=@(a,b,c,d,e) ((a*b)-(c*d))/e;
B_val=@(a,b,c,d,e) ((a*b)-(c*d))/e;
dlta_val=@(a,b,c) ((a*b)-(c^2));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part1 -  Best fit  line in the form of (y= b0 + (b1*t))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workout the bestfit line formula
AB_val = polyfit(t,y,2);
B0=AB_val(3);
B1=AB_val(2);

% Display the formula
text =sprintf('Part1:\nThe bestfit formula is "y = %g + (%g*(t))"\nB0 = %g\nB1 = %g\n\n',B0,B1,B0,B1);
disp(text);

% Dump the data to a txt file
textinput = sprintf('%s', text);
dlmwrite('output_uncertainty.txt',textinput,'delimiter',...
        '','roffset',1);

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part2 -  The uncertainty in measurement of 'y' (sigma_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error in 'y' (sigma_y)
for i = 1:c
    yi(i)= (y(i)-AB_val(3)-AB_val(2)*t(i))^2;
end
sigma_y = sqrt((1/(c-2))*sum(yi));

% Display the error in y (sigma_y)
text =sprintf('Part2:\nThe uncertainty in y =  %g\n\n',sigma_y);
disp(text);

% Dump the data to the file output_uncertainty.txt
textinput = sprintf('%s', text);
dlmwrite('output_uncertainty.txt',textinput,'-append','delimiter',...
        '','roffset',1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part3 -  The uncertainty in co-effeciency B0 and B1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delta
sgm_t1=0; sgm_t2=0;
for i = 1:c
    sgm_t1 = sgm_t1 + t(i)^2;
    sgm_t2 = sgm_t2 + t(i);
end
dlta= (c*sgm_t1)- ((sgm_t2)^2);
% Uncentainty in (A)
sigma_b0=sigma_y* sqrt(sgm_t1/dlta);
% Uncentainty in (B)
sigma_b1= sigma_y* sqrt(c/dlta);

% Display the uncertainty in A and B
text =sprintf('Part3:\nThe uncertainty in B0 =  %g\nThe uncertainty in B1 =  %g\n\n',sigma_b0,sigma_b1);
disp(text);

% Dump the data to the file output_uncertainty.txt
textinput = sprintf('%s', text);
dlmwrite('output_uncertainty.txt',textinput,'-append','delimiter',...
        '','roffset',1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part4 -  Find the best quadratic polynomial (y=B2+B3t+B4t^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amatrix = zeros(3,3);
Amatrix = [ (t.*t)' t' ones(length(t),1)];
Values = (((transpose(Amatrix)*Amatrix))^-1)*transpose(Amatrix)*y';

% Display the quadratic formula form
text =sprintf('Part4:\nThe quadratic formula for the best fit is "y =  (%g) + (%gt) + (%gt^2)"',...
    Values(3),Values(2),Values(1));
disp(text);

% Dump the data to a text file 
textinput = sprintf('%s', text);
dlmwrite('output_uncertainty.txt',textinput,'-append','delimiter',...
        '');
% Display the values
text =sprintf('B2 = %g\nB3 = %g\nB4 = %g\n\n',...
    Values(3),Values(2),Values(1));
disp(text);

% Dump the data to the file 
textinput = sprintf('%s', text);
dlmwrite('output_uncertainty.txt',textinput,'-append','delimiter',...
        '','roffset',1);
    
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part5 -  The uncertainty of y with a quadratic equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error in 'y' (sigma_y)
yi=0;
for i = 1:c
    yi(i)= (y(i)-Values(3)-Values(2)*t(i)-Values(1)*t(i)^2)^2;
end
sigma_y2 = sqrt((1/(c-2))*sum(yi));

% Display the error in y (sigma_y)
text =sprintf('Part5:\nThe uncertainty in y for a quadratic fit =  %g\n\n',sigma_y2);
disp(text);

% Dump the data to a txt file
textinput = sprintf('%s', text);
dlmwrite('output_uncertainty.txt',textinput,'-append','delimiter',...
        '','roffset',1);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part6 -  The uncertainty of the b3,b4,b5 values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build W matrix from previous uncertainty
errory = [(sigma_y2),(sigma_y2),(sigma_y2)];
inverseErrorSq = 1./(errory.*errory);
Wmatrix = diag(inverseErrorSq);
% Build Qmatrix
Q = (Amatrix'*Wmatrix*Amatrix)^(-1);
Q = sqrt(Q);

sigma_b2= Q(1);
sigma_b3= Q(5);
sigma_b4= Q(9);
% Display the uncertainty in b3,4,5
text =sprintf('Part6:\nThe uncertainty in B2 =  +/- %g\nThe uncertainty in B3 =  +/- %g\nThe uncertainty in B4 =  +/- %g\n\n',...
    sigma_b2,sigma_b3,sigma_b4);
disp(text);

% Dump the data to the file output_uncertainty.txt
textinput = sprintf('%s', text);
dlmwrite('output_uncertainty.txt',textinput,'-append','delimiter',...
        '','roffset',1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part7 -  Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the line fit plot onto the same graph
plot(t,y,'o');

% the line best fit
t1=linspace (t(1),t(end),1000);
bestfit=@(j) (AB_val(3)+(AB_val(2)*j));
plot(t1,bestfit(t1));

% the Quad of best fit
bestfit2=@(h) (Values(3)+(Values(2)*h)+ (Values(1)*h.^2));
plot(t1,bestfit2(t1));

legend('show')
legend('Original Data','Bestfit Data line','Bestfit Quadratic')
hold off;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part8 -  Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%workout the y value at t=5 and add them to the vectors
y(4)=bestfit(5);
t(4)= 5;

% REdo Error in 'y' (sigma_y) for a straight line
yi=0;
for i = 1:4
    yi(i)= (y(i)-AB_val(3)-AB_val(2)*t(i))^2;
end
sigma_y = sqrt((1/(4-2))*sum(yi));

% REdo Error in 'y' (sigma_y) for a straight line
y(4)=bestfit2(5);
yi=0;
for i = 1:4
    yi(i)= (y(i)-Values(3)-Values(2)*t(i)-Values(1)*t(i)^2)^2;
end
sigma_y2 = sqrt((1/(4-2))*sum(yi));