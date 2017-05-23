fid = fopen('asteroideja.txt', 'w');
tlist=dir('Dimensiot\');
tied=cell(length(tlist)-2, 1);
for i=3:length(tlist)
    tied{i-2}=strcat('Dimensiot\', tlist(i).name);
    fprintf(fid, tied{i-2});
end
fclose(fid);