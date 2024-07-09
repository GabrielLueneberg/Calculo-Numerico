x = randn(1000, 1);
fid = fopen('binario.bin', 'wb');
fwrite(fid, x,'double');
fclose(fid);

fid = fopen('ascii.txt', 'w');
numbytesascii = fprintf(fid,'%f', x);
fclose(fid);


%Quando o arquivo binario.bin é gerado, nota-se que seu tamanho é de 8kb.
%Isso está de acordo com o esperado visto que cada double precisa de 8bytes para ser representado, totalizando 8000 bytes.
%O arquivo txt tem tendencia a ser maior por cada caracter necessitar de 1 byte para ser representado
