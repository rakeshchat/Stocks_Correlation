if Input1==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data1 = importdata('Nikkei_165.csv');
    %Combine the date and time data together
    date_string = strcat(data1.textdata(2:end,1), {' '}, data1.textdata(2:end,2));
    %Create an array of time and date numbers for the X axis
    xdate = datenum(date_string, 'dd/mm/yy');
end
if Input1==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data1 = importdata('SP_194.csv');
    %Combine the date and time data together
    date_string = strcat(data1.textdata(2:end,1), {' '}, data1.textdata(2:end,2));
    %Create an array of time and date numbers for the X axis
    xdate = datenum(date_string, 'dd/mm/yy');
end