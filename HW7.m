%HW7
% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).  
% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.
% When X = 0 or 1; dX/dt = 0; X is not changing therefore at steady state. 

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 
a = 2;
syms X;
fplot(a*X*(1-X));

b = -2;
syms X;
fplot(b*X*(1-X));

% yes it depends on the sign of a. when a is positive, X = 0 is unstable X= 1 is stable; 
% when a is negative, the stability changes. 

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 
[timecourse1 nnmax] = timeSeries1(0.1, 2) ;
[timecourse1 nnmax] = timeSeries1(0.5, 3) ;


% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 
hold on 
for a = 0.1:0.1:4
    j = round(a*10)
Xt= rand(1,200);
for i =1:5000
    Xt_1 = a.*Xt.*(1-Xt);
    Xt = Xt_1; 
end 
%Xf(j,:) = Xt;

scatter(a*ones(1,200), Xt, 'r.');
end 
hold off 
% There is bistability at a = 3
% when a is small the system moves to zero
% when a is intermediate value system goes to a non-zero steady state. 
% when a is bigger than 3 system can go to many different points 


% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 

dA/dt = V/1+B^4 - A;
dB/dt = V/1+A^4 - B;

%
% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 
V = 5; 
rhs = @(t,x) [V/(1+x(2)^4)-x(1); V/(1+x(1)^4)-x(2)];
A01 = 2
B01 = 1

solx = ode23(rhs, [0 25],[A01, B01]);
figure (1)
plot(solx.x, solx.y(1,:), 'r.-'); hold on;
plot(solx.x, solx.y(2,:), 'g.-'); 
legend({'A','B'}) 

A02 = 1
B02 = 2
solx = ode23(rhs, [0 25],[A02, B02]);
figure (2)
plot(solx.x, solx.y(1,:), 'r.-'); hold on;
plot(solx.x, solx.y(2,:), 'g.-'); 
legend({'A','B'}) 

%
% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter. 
figure; hold on; 
 
for V = 0:0.05:5
    polycoeff = [1, -1,0,0,(1+V^4),-V] 
    rts = roots(polycoeff); 
    rts = rts(imag(rts)==0); 
    plot(V*ones(length(rts),1),rts,'r.');
end 

title('n =1')