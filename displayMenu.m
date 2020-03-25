function choice = displayMenu(options)
% DISPLAYMENU Displays a menu of options, ask the user to choose an item
% and returns the number of the menu item chosen.
%
% Usage: choice = displayMenu(options)
%
% Input    options   Menu options (cell array of strings)
% Output   choice    Chosen option (integer)
%
% Author: Mikkel N. Schmidt, mnsc@dtu.dk, 2015
% Display menu options
for i = 1:length(options)
    fprintf('%d. %s\n', i, options{i});
end
% Get a valid menu choice
choice = 0;
while ~any(choice == 1:length(options))
    choice = inputNumber('Please choose a menu item: ');
    if ~any(choice == 1:length(options))
        disp('Not valid')
    end
end
end