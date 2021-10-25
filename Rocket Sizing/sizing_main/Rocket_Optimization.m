%% Clearing Screen
clc
clear

%% Main

% Step 1- Request Inputs
[Vdesired, Vratio, Hmax, DHratioMax, payload] = reqInput();

% Step 2 - Initialize variables
dVlow = Vdesired * Vratio; % deltaV lower stage
dVup = Vdesired - dVlow; % deltaV upper stage
lengthL = linspace(0, Hmax, 100); % 100 entries for Height
ratioDL = linspace(0, DHratioMax, 1000); % 100 entries for Diameter/Height ratio
loopbroken = 0; % To create While loop
a = 0; % Indexing variable

% Step 3 - Rocket Optimization
while loopbroken == 0
    a = a + 1; % Indexing While Loop
    temp(a) = rand(); % Call First Stage
    temp2(a) = rand(); % Call Second stage
    temp3(a) = temp(a) + temp2(a);
    if temp3(a) < 0.1
        loopbroken = 1;
        fprintf("Done.\n")
        disp(temp3(a))
    else
        loopbroken = 0;
    end
end

%% Function Definitions
function [Vdesired, Vratio, Hmax, DHratioMax, payload] = reqInput()
prompt = 'What is the desired total deltaV? (km/s)  ';
Vdesired = input(prompt);

prompt = 'What is the desired deltaV ratio between two stages? (lower to upper)  ';
Vratio = input(prompt);

prompt = 'What is the height limit for either stage? (m)  ';
Hmax = input(prompt);

prompt = 'What is the max ratio for diameter / height?  ';
DHratioMax = input(prompt);

prompt = 'What is the initial payload?  ';
payload = input(prompt);
end

