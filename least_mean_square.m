%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:   
%        least mean square
%        error propigation page 196 (8.34) formula
%            sigma_A, sigma_B, Sigma_z 
%
% Author: Alison
% Created: Fall 2016
% Modfied: 6 Feb 2016
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
title ('Title of Some Sort','FontName','Arial','FontSize',20)
xlabel('time','FontName','Times','FontSize',18,'FontWeight','bold')
ylabel('z','FontName','Times','FontSize',18,'FontWeight','bold')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% values of x, and y we need to least mean square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = [10,40,70,100,130,160];
y = [188,102,60,18,16,5];
[~,c]=size(x);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formulas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%symbolics
syms v;

%formulas
A_val=@(a,b,c,d,e) ((a*b)-(c*d))/e;
B_val=@(a,b,c,d,e) ((a*b)-(c*d))/e;
dlta_val=@(a,b,c) ((a*b)-(c^2));
z=@(v) log(v);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve for error propigation page 196 (8.34) formula
sigma_z= diff(z,v) * sqrt(v);   %this shows the answer is 1/(v^(1/2))
text =sprintf('%s',sigma_z)

% Dump the data to a text file
textinput = sprintf('The uncertainty of (z) = %s', text);
dlmwrite('output_lms.txt',textinput,'delimiter',...
        '','roffset',1);

%workout the error for the error bars    
for i= 1:c
err(i)=(1/(sqrt(y(i))));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part b matrix weighted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z=0; 
for i= 1:c
z(i) = log(y(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%delta value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1=0; part2=0; part3=0; part4=0;
for i= 1:c
part1= part1 + y(i);
part2= part2 + (y(i)*(x(i)^2));
part3= part3 + (y(i)*x(i));
end
dlta= dlta_val(part1,part2,part3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (A) value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1=0; part2=0; part3=0; part4=0;
for i= 1:c
part1= part1 + (y(i)*(x(i)^2));
part2= part2 + ((y(i))*z(i));
part3= part3 + (y(i)*x(i));
part4= part4 + (y(i)*x(i)*z(i));
end
A= A_val(part1,part2,part3,part4,dlta);
%uncentainty in (A)
sigma_A=sqrt(part1/dlta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (B) value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1=0; part2=0; part3=0; part4=0;
for i= 1:c
part1= part1 + y(i);
part2= part2 + (z(i)*x(i)*y(i));
part3= part3 + (y(i)*x(i));
part4= part4 + (y(i)*z(i));
end
B=B_val(part1,part2,part3,part4,dlta);
%uncentainty in (A)
sigma_B=sqrt(part1/dlta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make table of t,v,z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text=sprintf('t-value\t\tV-value\t\tz-value');
disp(text);
textinput= sprintf('%s',text);
dlmwrite('output_lms.txt',textinput,'-append','delimiter',...
        '','roffset',1);

for i = 1:c
text=sprintf('%g\t\t\t%g\t\t\t%g',x(i),y(i),z(i));
disp(text);
textinput= sprintf('%s',text);
dlmwrite('output_lms.txt',textinput,'-append','delimiter',...
        '');
end
text=sprintf('\n\nSigma_A = %g',sigma_A);
disp(text);
textinput= sprintf('%s',text);
dlmwrite('output_lms.txt',textinput,'-append','delimiter',...
        '','roffset',1);
text=sprintf('Sigma_B = %g',sigma_B);
disp(text);
textinput= sprintf('%s',text);
dlmwrite('output_lms.txt',textinput,'-append','delimiter',...
        '','roffset',1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part C graph of z vs t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot z vs t

errorbar(x,z,err,'.')

% best fit values for z vs t
AB_val = polyfit(x,z,2);

%Formula of the line fit plot
f=@(c) (AB_val(3)+(AB_val(2)*c));
x1=linspace (x(1),x(end),1000);

%Plot the line fit plot onto the same graph
plot(x1,f(x1));
legend('show')
legend('Original Data','Bestfit Data')
hold off;