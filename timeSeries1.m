function [timecourse1 nnmax] = timeSeries1(x0, P) 
rhs = @(t,x) P*x*(1-x)
sol = ode23(rhs, [0 10],x0)
timecourse1 = sol.y ;
timecourse_round=round(timecourse1, 2) ;
k = find(timecourse_round >= 0.99);
nnmax = sol.x(k(1)) ;