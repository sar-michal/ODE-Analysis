clearvars
obj1_xy = [0.5, sqrt(3)/2];
obj2_xy = [0.5, -sqrt(3)/2];
obj3_xy = [-1, 0];
obj1_dt = [sqrt(3)/2, -1/2];
obj2_dt = [-sqrt(3)/2, -1/2];
obj3_dt = [0, 1];
z0 = [obj1_xy, obj2_xy, obj3_xy, obj1_dt, obj2_dt, obj3_dt];
tspan = [0, 2*pi];

figure(1)
clf
hold on;
yscale('log');
xscale('log');
h = logspace(-3,0,256);

%euler explicit
z_euler_delta = zeros(1, length(h));
for i = 1:length(h)
    t_steps = tspan(1):h(i):tspan(2);
    z_euler = zeros(length(z0), length(t_steps));
    z_euler(:,1) = z0;
    for n = 2:(length(t_steps))
        z_euler(:, n) = z_euler(:,n-1) + h(i)*(odeFun(t_steps(n-1), z_euler(:,n-1)));
    end
    z_euler_y2_ref = zeros(1, length(t_steps));
    z_euler_y2_ref(1, :) = cos(t_steps + (5/6)*pi);
    z_euler_delta(i) = 1/(length(t_steps))*sum((z_euler(4, :) - z_euler_y2_ref).^2);
end
plot(h, z_euler_delta);

%Adams-Bashforth of order 2
z_adams_delta = zeros(1, length(h));
for i = 1:length(h)
    t_steps = tspan(1):h(i):tspan(2);
    z_adams = zeros(length(z0), length(t_steps));
    z_adams(:,1) = z0;
    z_adams(:, 2) = z_adams(:,1) + h(i)*(odeFun(t_steps(1), z_adams(:,1)));
    for n = 3:(length(t_steps))
        z_adams(:, n) = z_adams(:,n-1) + (h(i)/2)*(3*odeFun(t_steps(n-1), z_adams(:,n-1)) - odeFun(t_steps(n-2), z_adams(:,n-2)));
    end
    z_adams_y2_ref = zeros(1, length(t_steps));
    z_adams_y2_ref(1, :) = cos(t_steps + (5/6)*pi);
    z_adams_delta(i) = 1/(length(t_steps))*sum((z_adams(4, :) - z_adams_y2_ref).^2);
end
plot(h, z_adams_delta);

%Gear explicit of order 2
z_gear_delta = zeros(1, length(h));
for i = 1:length(h)
    t_steps = tspan(1):h(i):tspan(2);
    z_gear = zeros(length(z0), length(t_steps));
    z_gear(:,1) = z0;
    z_gear(:, 2) = z_gear(:,1) + h(i)*(odeFun(t_steps(1), z_gear(:,1)));
    for n = 3:(length(t_steps))
        z_gear(:, n) = z_gear(:,n-2) + h(i)*(2*odeFun(t_steps(n-1), z_gear(:,n-1)));
    end
    z_gear_y2_ref = zeros(1, length(t_steps));
    z_gear_y2_ref(1, :) = cos(t_steps + (5/6)*pi);
    z_gear_delta(i) = 1/(length(t_steps))*sum((z_gear(4, :) - z_gear_y2_ref).^2);
end
plot(h, z_gear_delta);
ylabel('\Delta_{y_2}', 'FontName', 'Cambria Math');
xlabel('h', 'FontName', 'Cambria Math');
legend('Euler Explicit Method', 'Adams-Bashforth Method', 'Gear Explicit Method', 'FontName', 'Times New Roman', 'Location','southeast');

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

    a1 = (-sqrt(3))*(xy1 - xy2)/(r12)^3 - (sqrt(3))*(xy1 - xy3)/(r31)^3;
    a2 = (-sqrt(3))*(xy2 - xy3)/(r23)^3 - (sqrt(3))*(xy2 - xy1)/(r12)^3;
    a3 = (-sqrt(3))*(xy3 - xy1)/(r31)^3 - (sqrt(3))*(xy3 - xy2)/(r23)^3;

    vector = [dt1; dt2; dt3; a1; a2; a3];
end