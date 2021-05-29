clear;
format long;
critical_temperature = 631.15; %k
critical_pressure = 3.2; %Mpa
critical_volume = 0.4370;
omega = 0.3377; 
T = [400:5:620];
obj = srk(critical_temperature, critical_pressure, critical_volume , omega);
obj = obj.calculate_constant_for_EOS();

psat = [];
vol_l = [];
vol_g = [];
length = size(T);
for i = 1:length(2)
    obj = obj.calculate_alpha(T(i));
    %obj.plot(T(i))
    psat = [psat obj.find_psat(T(i))];
    volume = sort(obj.find_volume(psat(i) , T(i)));
    vol_l = [vol_l volume(1)];
    vol_g = [vol_g volume(3)];
end
[pc,tc,vc,upper_dome_region]= obj.find_critical_temp_and_pressure();
fprintf ("critical pressure obtained from the equation is %.4f\n",pc);
fprintf ("critical temperature obtained from the equation is %.4f\n",tc);
fprintf ("critical volume obtained from the equation is %.4f\n",vc);
hold on;
%graph is only plotted for molar volume 0 to 9 litres
%plotting lower dome region 
plot(vol_l ,psat,'b--o')
s = size(vol_g);
plot(vol_g(:,20:s(2)) ,psat(:,20:s(2)),'b--o')
ylabel("Pressure (Mpa)");
xlabel("molar volume (L/mol)");
title("PV CURVE OF SRK FOR CUMENE ");
%fitting a cubic curve in upper_dome_region
upper_dome_vol = [upper_dome_region(:,3) ;upper_dome_region(:,4)];
upper_dome_psat = [upper_dome_region(:,2) ;upper_dome_region(:,2)];
p = polyfit(upper_dome_vol,upper_dome_psat,3);
%plotting upper dome region
x = linspace(vol_l(s(2)),vol_g(s(2))-0.15, 20);
y = polyval(p,x);
x = [vol_l(s(2)) x vol_g(s(2))];
y = [psat(s(2)) y psat(s(2))];
plot(x,y,'m--')
x = linspace(vol_l(s(2)),vol_g(s(2))-0.15,5);
y = polyval(p,x);
x = [vol_l(s(2)) x vol_g(s(2))];
y = [psat(s(2)) y psat(s(2))];
plot(x,y,'mo')
%plotting critical point on the curve obtained from the literature
plot(critical_volume,critical_pressure,'c*','LineWidth', 1)
%plotting critical point on the curve obtained by solving the equation of state 
plot(vc,pc,'k*','LineWidth', 1)
%plotting curves in supercritical region
T = [critical_temperature:10:critical_temperature+30];
length = size(T);
for i = 1 :length(2)
    obj = obj.calculate_alpha(T(i));
    obj.plot(T(i),obj.b,9);
end
%plotting curves in subcritical region
T = [520 540 560 580 600 620];

for i = 1:6
    obj = obj.calculate_alpha(T(i));
    psat = obj.find_psat(T(i));
    volume = sort(obj.find_volume(psat , T(i)));
    vol_l = volume(1);
    vol_g = volume(3);
    obj.plot(T(i),obj.b,vol_l);
    plot([vol_l vol_g], [psat psat],'r','LineWidth',1);
    obj.plot(T(i),vol_g,9);
end
