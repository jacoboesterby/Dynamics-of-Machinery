function num = inputNumber(prompt)
% INPUTNUMBER Prompts user to input a number
%
% Usage: num = inputNumber(prompt) Displays prompt and asks user to input a
% number. Repeats until user inputs a valid number.
%
% Author: Mikkel N. Schmidt, mnsc@dtu.dk, 2015

while true    
    num = str2double(input(prompt, 's'));
    if ~isnan(num)
        break;
    end
end
