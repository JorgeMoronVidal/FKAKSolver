%
% utility to read MLMC data from input text file
%
% mlmc_read(filename)
%

function [del1 del2 var1 var2 kur1 chk1 cost  ...
          Eps mlmc_cost std_cost ls Nls ]  =  mlmc_read(filename)

%
% read in data
%

fid = fopen([filename '.txt'],'r');

line = '    ';
while (length(line)<20) | (strcmp(line(1),'-')==0)
  line = [ fgetl(fid) '    ' ];
end

line = fgetl(fid);
l    = 1;
while (length(line)>10)
  data = sscanf(line,'%f');
  del1(l) = data(2);
  del2(l) = data(3);
  var1(l) = data(4);
  var2(l) = data(5);
  kur1(l) = data(6);
  chk1(l) = data(7);
  cost(l) = data(8);

  line = fgetl(fid);
  l    = l+1;
end

L = l-2;

line = '    ';
while (length(line)<20) | (strcmp(line(1),'-')==0)
  line = [ fgetl(fid) '    ' ];
end

line = fgetl(fid);
l    = 1;

while (length(line)>10)
  data = sscanf(line,'%f');
  Eps(l)       = data(1);
  mlmc_cost(l) = data(3);
  std_cost(l)  = data(4);
  len          = length(data)-5;
  ls(1:len,l)  = 0:len-1;
  Nls(1:len,l) = data(6:end);

  line = fgetl(fid);
  l    = l+1;
end

fclose(fid);

