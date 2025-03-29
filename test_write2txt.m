% Define the file name
fileName = 'test_results/example.txt';

% Open the file for writing
fileID = fopen(fileName, 'w');

% Check if the file opened successfully
if fileID == -1
    error('Cannot open file for writing.');
end

% Define the lines to write
lines = {
    'This is the first line.';
    'This is the second line.';
    'This is the third line.'
};

% Write each line to the file
for i = 1:length(lines)
    fprintf(fileID, '%s hell0o %d\n', lines{i}, 1);
end
x = 1:10;


% Close the file
fclose(fileID);

disp(['File written successfully to ', fileName]);
