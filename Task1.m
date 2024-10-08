clearvars
%{
m1 = 1;
m2 = 1;
m3 = 1;
G = 1;
%}
obj1_xy = [0.0132604844, 0];
obj2_xy = [1.4157286016, 0];
obj3_xy = [-1.4289890859, 0];
obj1_dt = [0, 1.0541519210];
obj2_dt = [0, -0.2101466639];
obj3_dt = [0, -0.8440052572];
z0 = [obj1_xy, obj2_xy, obj3_xy, obj1_dt, obj2_dt, obj3_dt];
tspan = [0, 32.584945];

%ode45
options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[t, z] = ode45(@odeFun, tspan, z0, options);
figure(1)
clf
hold on;
plot(z(:,1), z(:, 2));
plot(z(:,3), z(:, 4));
plot(z(:,5), z(:, 6));
ylabel('y', 'FontName', 'Cambria Math');
xlabel('x', 'FontName', 'Cambria Math');
legend('Object 1', 'Object 2', 'Object 3', 'FontName', 'Times New Roman');
axis equal
hold off;

%euler explicit
h = 0.001;
t_steps = tspan(1):h:tspan(2);
z_euler = zeros(length(z0), length(t_steps));
z_euler(:,1) = z0;
for n = 2:(length(t_steps))
    z_euler(:, n) = z_euler(:,n-1) + h*(odeFun(t_steps(n-1), z_euler(:,n-1)));
end
figure(2)
clf
hold on;
plot(z_euler(1, :), z_euler(2, :));
plot(z_euler(3, :), z_euler(4, :));
plot(z_euler(5, :), z_euler(6, :));
ylabel('y', 'FontName', 'Cambria Math');
xlabel('x', 'FontName', 'Cambria Math');
legend('Object 1', 'Object 2', 'Object 3', 'FontName', 'Times New Roman');
axis equal
hold off;

%Adams-Bashforth of order 2
z_adams = zeros(length(z0), length(t_steps));
z_adams(:,1) = z0;
z_adams(:, 2) = z_adams(:,1) + h*(odeFun(t_steps(1), z_adams(:,1)));
for n = 3:(length(t_steps))
    z_adams(:, n) = z_adams(:,n-1) + (h/2)*(3*odeFun(t_steps(n-1), z_adams(:,n-1)) - odeFun(t_steps(n-2), z_adams(:,n-2)));
end
figure(3)
clf
hold on;
plot(z_adams(1, :), z_adams(2, :));
plot(z_adams(3, :), z_adams(4, :));
plot(z_adams(5, :), z_adams(6, :));
ylabel('y', 'FontName', 'Cambria Math');
xlabel('x', 'FontName', 'Cambria Math');
legend('Object 1', 'Object 2', 'Object 3', 'FontName', 'Times New Roman');
axis equal
hold off;

%Gear explicit of order 2
z_gear = zeros(length(z0), length(t_steps));
z_gear(:,1) = z0;
z_gear(:, 2) = z_gear(:,1) + h*(odeFun(t_steps(1), z_gear(:,1)));
for n = 3:(length(t_steps))
    z_gear(:, n) = z_gear(:,n-2) + h*(2*odeFun(t_steps(n-1), z_gear(:,n-1)));
end
figure(4)
clf
hold on;
plot(z_gear(1, :), z_gear(2, :));
plot(z_gear(3, :), z_gear(4, :));
plot(z_gear(5, :), z_gear(6, :));
ylabel('y', 'FontName', 'Cambria Math');
xlabel('x', 'FontName', 'Cambria Math');
legend('Object 1', 'Object 2', 'Object 3', 'FontName', 'Times New Roman');
axis equal
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

    a1 = (-1)*(xy1 - xy2)/(r12)^3 - (xy1 - xy3)/(r31)^3;
    a2 = (-1)*(xy2 - xy3)/(r23)^3 - (xy2 - xy1)/(r12)^3;
    a3 = (-1)*(xy3 - xy1)/(r31)^3 - (xy3 - xy2)/(r23)^3;

    vector = [dt1; dt2; dt3; a1; a2; a3];
end