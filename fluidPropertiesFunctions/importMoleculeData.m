function [const_up,const_low] = importMoleculeData(gasMolecule)
% Function to extract NASA constants for polynomial description of cp,
% gamma, R
% INPUT
% gasMolecule: string with name of the molecule
%
% OUTPUT
% const_up,const_low: vectors with NASA constants of the type:
% [a0 a1 a2 a3 a4]

dati='NASA_thermoData_1999.txt';
text = fileread(dati);
textLines = splitlines(text);
expression=horzcat('^',gasMolecule,'\s');
searchOutput=(regexp(textLines,expression));

% Convert cell to vector i replace empty with Nan
searchOutput(cellfun('isempty',searchOutput)) = {NaN};  
% Vectorial output
searchOutput = [searchOutput{:}];   
lineData=find(searchOutput==1);

%Data for the desired molecule
DataLines=textLines(lineData+1:lineData+3);           
line1=(mat2str(cell2mat(DataLines(1))));
line2=(mat2str(cell2mat(DataLines(2))));
line3=(mat2str(cell2mat(DataLines(3))));

const_up=[str2double(line1(2:16)) str2double(line1(17:31)) str2double(line1(32:46)) str2double(line1(47:61)) str2double(line1(62:76)) str2double(line2(2:16)) str2double(line2(17:31))];
const_low=[str2double(line2(32:46)) str2double(line2(47:61)) str2double(line2(62:76)) str2double(line3(2:16)) str2double(line3(17:31)) str2double(line3(32:46)) str2double(line3(47:61))];

end

