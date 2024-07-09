x = 0;
tolerancia = 0.0000001;

for i = 1:9
    x = x + 1/3;
     fprintf('o valor de x: %f\n', x)
    if abs(x - 1) < tolerancia
        disp('chegamos em x = 1');
    elseif abs(x - 2) < tolerancia
        disp('chegamos em x = 2');
    elseif abs(x - 3) < tolerancia
        disp('chegamos em x = 3');
    end
end
