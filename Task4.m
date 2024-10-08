clearvars;
table = readtable("data_42.csv");
obj1_xy = [table.x1(1), table.y1(1)];
obj2_xy = [table.x2(1), table.y2(1)];
obj3_xy = [table.x3(1), table.y3(1)];
t_change = table.t(2)-table.t(1);
obj1_dt = [(table.x1(2)-table.x1(1))/t_change, (table.y1(2)-table.y1(1))/t_change];
obj2_dt = [(table.x2(2)-table.x2(1))/t_change, (table.y2(2)-table.y2(1))/t_change];
obj3_dt = [(table.x3(2)-table.x3(1))/t_change, (table.y3(2)-table.y3(1))/t_change];
p_ini = [obj1_dt, obj2_dt, obj3_dt];
z0 = [obj1_xy, obj2_xy, obj3_xy, p_ini];

options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[t, z] = ode45(@odeFun, table.t, z0, options);

%parameters = fminsearch(@minimisation, p_ini);
%Approximated parameters
parameters = [0.725359552559566,	0.00343891969394009,	-0.373901025070656,	-0.629540103929511,	-0.371860831386569,	0.640514560568173];
z0_approx = [obj1_xy, obj2_xy, obj3_xy, parameters];
[t, z_approx] = ode45(@odeFun, table.t, z0_approx, options);

figure(1);
clf;
hold on;
axis equal
p = plot(table.x1, table.y1, '.');
p.MarkerSize = 8;
plot(z(:,1), z(:, 2), LineWidth = 1.5);
plot(z_approx(:,1), z_approx(:, 2), LineWidth = 1.5);
ylabel('y', 'FontName', 'Cambria Math');
xlabel('x', 'FontName', 'Cambria Math');
legend('Given coordinates', 'Simple estimation', 'fminsearch estimation', 'FontName', 'Times New Roman', 'Location', 'northwest');


hold off;

figure(2);
clf;
hold on;
axis equal
p = plot(table.x2, table.y2, '.');
p.MarkerSize = 8;
plot(z(:,3), z(:, 4), LineWidth = 1.5);
plot(z_approx(:,3), z_approx(:, 4), LineWidth = 1.5);
ylabel('y', 'FontName', 'Cambria Math');
xlabel('x', 'FontName', 'Cambria Math');
legend('Given coordinates', 'Simple estimation', 'fminsearch estimation', 'FontName', 'Times New Roman', 'Location', 'northwest');
hold off;

figure(3);
clf;
hold on;
axis equal
p = plot(table.x3, table.y3, '.');
p.MarkerSize = 8;
plot(z(:,5), z(:, 6), LineWidth = 1.5);
plot(z_approx(:,5), z_approx(:, 6), LineWidth = 1.5);
ylabel('y', 'FontName', 'Cambria Math');
xlabel('x', 'FontName', 'Cambria Math');
legend('Given coordinates', 'Simple estimation', 'fminsearch estimation', 'FontName', 'Times New Roman', 'Location', 'northwest');
hold off;

function vector = odeFun(t, w) 
    xy1 = w(1:2);
    xy2 = w(3:4);
    xy3 = w(5:6);
    dt1 = w(7:8);
    dt2 = w(9:10);
    dt3 = w(11:12);
    r12=norm(xy1 - xy2);
    r31=norm(xy3 - xy1);
    r23=norm(xy2 - xy3);

    a1 = (-1)*(xy1 - xy2)/(r12)^3 - (1)*(xy1 - xy3)/(r31)^3;
    a2 = (-1)*(xy2 - xy3)/(r23)^3 - (1)*(xy2 - xy1)/(r12)^3;
    a3 = (-1)*(xy3 - xy1)/(r31)^3 - (1)*(xy3 - xy2)/(r23)^3;

    vector = [dt1; dt2; dt3; a1; a2; a3];
end

function J = minimisation(p)
    table = readtable("data_42.csv");
    xy1_ini = [table.x1(1), table.y1(1)];
    xy2_ini = [table.x2(1), table.y2(1)];
    xy3_ini = [table.x3(1), table.y3(1)];
    initial_parameters = [xy1_ini, xy2_ini, xy3_ini, p];
    options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
    [t, z] = ode45(@odeFun, table.t, initial_parameters, options);
    xy1 = z(:, 1:2);
    xy2 = z(:, 3:4);
    xy3 = z(:, 5:6);
    
    J = sum(norm(xy1 - [table.x1, table.y1])) + sum(norm(xy2 - [table.x2, table.y2])) + sum(norm(xy3 - [table.x3, table.y3]));
        
end