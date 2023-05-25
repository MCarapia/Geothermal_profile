function[T_sol, T_liq] = T_sol_liq (Profundidad_km)

MP=Profundidad_km.*27; % Conversion de kilometros a Megapascales
GP=MP./1000; % Conversion de Megapascales a Gigapascales

for i=1:length(GP)
    if GP(i)<=2.7
        Tsol(:,i) = 1120.661 + 273.15 + (132.899*GP(i)) - (5.904*(GP(i)^2));
    elseif (2.7 < GP(i)) && (GP(i)<= 22.5)
        Tsol(:,i) = 1086 + 273.15  - 5.7*GP(i) + 390*log(GP(i));
    elseif GP(i) > 22.5
        Tsol(:,i) = 1762.722 + 31.595*GP(i) - 0.102*GP(i)^2;
    else
    end
end


for i=1:length(GP)
    if GP(i) < 22.5
        Tliq(:,i) =  2014.497 + (37.743*GP(i)) - (0.472*(GP(i)^2));
    elseif GP(i) >= 22.5
        Tliq(:,i) = 1803.547 + 50.810*GP(i) - 0.185*GP(i)^2;
    else
    end
end


T_sol = Tsol;
T_liq = Tliq; 
Profundidad = 6370- Profundidad_km;


end


