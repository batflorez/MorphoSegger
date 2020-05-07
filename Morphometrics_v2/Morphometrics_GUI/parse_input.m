%Tristan Ursell
%Oct 2013
%
%
% text_in='XXX, YYY,ZZZ  ';
% name_struct{1}='XXX';
% name_struct{2}='YYY';
% name_struct{3}='ZZZ';
%
% default is to parse names separated by commas
%
% Can also specify a parsing character by:
%
% name_struct=parse_input(text_in,'*');
%

function name_struct=parse_input(text_in,varargin)
if isempty(varargin)
    parse_char=',';
else
    parse_char=varargin{1};
end

text_in=strtrim(text_in);
%parse channel name inputs
if isempty(text_in)
    name_struct={};
else
    remain=text_in;
    i=1;
    
    while ~isempty(remain)
        [token, remain] = strtok(remain, parse_char);       
        name_struct{i}=strtrim(token);       
        i=i+1;
    end
end