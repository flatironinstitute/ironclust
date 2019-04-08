n_e = 8;

e = [sym('e', [2,n_e]); zeros(1,n_e)]; % electrode position
x = sym('x',[3,1]);
v = sym('v',[1,n_e]);
v_fit = 1./sqrt(sum((e-x).^2));
obj = 1 - v*v_fit/abs(v)/abs(v_fit);

simplify(diff(obj, a))
simplify(diff(obj, x(1)))
simplify(diff(obj, x(2)))
simplify(diff(obj, x(3)))

%%


(y1 - a*(x3^2 + (e1_1 - x1)^2 + (e2_1 - x2)^2))*(x3^2 + (e1_1 - x1)^2 + (e2_1 - x2)^2) - 2*(y2 - a*(x3^2 + (e1_2 - x1)^2 + (e2_2 - x2)^2)) * (x3^2 + (e1_2 - x1)^2 + (e2_2 - x2)^2) - 2*(y3 - a*(x3^2 + (e1_3 - x1)^2 + (e2_3 - x2)^2))*(x3^2 + (e1_3 - x1)^2 + (e2_3 - x2)^2) - 2*(y4 - a*(x3^2 + (e1_4 - x1)^2 + (e2_4 - x2)^2))*(x3^2 + (e1_4 - x1)^2 + (e2_4 - x2)^2)
(y1 - a*(x3^2 + (e1_1 - x1)^2 + (e2_1 - x2)^2))*(2*e1_1 - 2*x1) + 2*a*(y2 - a*(x3^2 + (e1_2 - x1)^2 + (e2_2 - x2)^2))*(2*e1_2 - 2*x1) + 2*a*(y3 - a*(x3^2 + (e1_3 - x1)^2 + (e2_3 - x2)^2))*(2*e1_3 - 2*x1) + 2*a*(y4 - a*(x3^2 + (e1_4 - x1)^2 + (e2_4 - x2)^2))*(2*e1_4 - 2*x1)
(y1 - a*(x3^2 + (e1_1 - x1)^2 + (e2_1 - x2)^2))*(2*e2_1 - 2*x2) + 2*a*(y2 - a*(x3^2 + (e1_2 - x1)^2 + (e2_2 - x2)^2))*(2*e2_2 - 2*x2) + 2*a*(y3 - a*(x3^2 + (e1_3 - x1)^2 + (e2_3 - x2)^2))*(2*e2_3 - 2*x2) + 2*a*(y4 - a*(x3^2 + (e1_4 - x1)^2 + (e2_4 - x2)^2))*(2*e2_4 - 2*x2)
(y1 - a*(x3^2 + (e1_1 - x1)^2 + (e2_1 - x2)^2)) - 4*a*x3*(y2 - a*(x3^2 + (e1_2 - x1)^2 + (e2_2 - x2)^2)) - 4*a*x3*(y3 - a*(x3^2 + (e1_3 - x1)^2 + (e2_3 - x2)^2)) - 4*a*x3*(y4 - a*(x3^2 + (e1_4 - x1)^2 + (e2_4 - x2)^2))
 

%%
obj = sum((v-v_fit).^2);
solve([diff(obj,x(1)); diff(obj,x(2)); diff(obj,x(3)); diff(obj,a)] == [0;0;0;0], {x,a})


%%
prob = optimproblem('ObjectiveSense','max');
x = optimvar('x',2,1,'LowerBound',0);
prob.Objective = x(1) + 2*x(2);
