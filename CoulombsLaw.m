% INITIALIZE MATLAB
clc;
clear;
close all;


% CHARGE SPACE containing all the Charges
global charge_space; global charge_space_permittivity;
charge_space_permittivity = 8.85418782 * 10^(-12);
charge_space = struct('x_coord', 'N', 'y_coord', [], 'z_coord', [], 'mag', [], 'col', 'g');


% BACK END FUNCTIONS
% Charge Space Manipulation
clear_charge_space = @clear_space;  % Removes all charges from Charge Space
charge_place = @charge_create;  % Creates a charge at coordinates pointed by parameters 'x, y, z' in charge space of magnitude given by parameter 'mag' and plots it if parameter 'show' is 1
display_charge = @display_charge_properties;    % Display location and magnitude of charge given by parameter 'obj' in charge space
charge_destroy = @charge_dest;  % Removes charge given by parameter 'obj' from charge space
n_random_charges = @create_rnd_charges; % Generate number of random charges as given by parameter 'number'
% Locating Charges
get_chargenum_bycoord = @find_charge_bycoord;   % Get number of Charge pointed by given coordinate parameters 'x, y, z' in Charge Space
get_charge_bynum = @find_charge;       % Get Charge object pointed by parameter 'charge_number' in Charge Space 
% Finding physical Force and Electric Field Strength values
charge_pair_force = @force_on_two;  % Calculates the three dimensional components of force on charge parameter 'obj' by charge parameter 'obj2'
charge_net_force = @net_force_on;   % Superposition Theorem: Calculates the three dimensionalc components of the net force on charge parameter 'obj' by all other Charges
net_field_at = @net_field_on;       % Superposition Theorem: Calculates the three dimensionalc components of the net Electric Field strength at coordinate parameters 'x, y, z' in Charge Space

% CLIENT FUNCTIONS
plot_space = @plot_ch;  % Function that creates the Charge Space Menu
menu = @help;   % Function that creates the Calculate Values with Coulomb's Law menu
calc = @calculator; % Function that creates the Calculate Values with Coulomb's Law menu
plot_graphs = @custom_grapher;


% --------------------------- MAIN STARTS HERE -------------------------------
disp("This is a MATLAB program to simulate Electrostatic phenomenon, particularly Coulomb's Law.");
disp("We maintain a simulated 'Charge space' assumed to be an infinite homogenous insulating medium.");
disp("The assumption is that every charged body in it remains stationary.")
fprintf('\n');
info_text();
%  --------------------------- MAIN ENDS HERE -------------------------------


% FUNCTION DEFINITIONS
% Function that creates the Charge Space Menu
function help()
    global charge_space; choice = 0;
    while choice ~= 6
        clc;
        fprintf("\t\t\t\t\t\t\t\t\t\tCHARGE SPACE MENU");
        if (charge_space(1).x_coord == 'N')
            n_charges = 0;
        else
            n_charges = length(charge_space);
        end
        fprintf("\n")
        disp(['The Charge space currently has ', num2str(n_charges), ' charges. Exit menu to interact with Charge space plot.']);
        disp("Enter 1 to place a specific new charge into Charge Space: ");
        disp("Enter 2 to remove specific charge from Charge Space or 7 to remove all charges: ");
        disp("Enter 3 to generate N random charges in Charge Space: ");
        disp("Enter 4 to display properties of a specific charge in Charge Space: ");
        disp("Enter 5 to exit Menu and show 3 Dimensional plot of current Charge Space: ");
        choice = input("==> ");
        clc;
        switch choice
            case 1
                disp("Enter the x, y, z coordinates of the charge and it's magnitude: ");
                x = input("X coordinate: ");
                y = input("Y coordinate: ");
                z = input("Z coordinate: ");
                mag = input("Magnitude: ");
                charge_create(x, y, z, mag);
                input("Charge placed. Press ENTER to continue!");

            case 2
                ch = 0;
                while ch ~= 1 && ch ~= 2
                    ch = input("Enter 1 to remove charge by charge number, 2 to remove charge by it's coordinate: ");                   
                    fprintf('\n');
                    if ch == 1
                        rem_num = input("Enter the charge number to remove: ");
                        charge_dest(rem_num);
                        break;
                    elseif ch == 2
                        disp("Enter the x, y, z coordinates of the charge:"); 
                        x = input("X coordinate: ");
                        y = input("Y coordinate: ");
                        z = input("Z coordinate: ");
                        charge_dest(find_charge_bycoord(x, y, z));
                        break;
                    end
                end

            case 3
                gen_num = input("Enter number of random charges to generate in Charge Space: ");
                create_rnd_charges(gen_num);
                fprintf('\n');
                disp([num2str(gen_num), ' charges generated! Press Enter to continue.']);
                input("");

            case 4
                ch = 0;
                while ch ~= 1 && ch ~= 2
                    disp("Enter 1 to show information of Charge pointed by it's number in Charge Space: ");
                    disp("Enter 2 to show information of Charge pointed by it's coordinate in Charge Space: ");
                    ch = input("===>");
                    if ch == 1
                        rem_num = input("Enter the charge's number whose information is to be viewed: ");
                        display_charge_properties(find_charge(rem_num));
                        input('');
                    elseif ch == 2
                        disp("Enter the x, y, z coordinates of the charge:");
                        x = input("X coordinate: ");
                        y = input("Y coordinate: ");
                        z = input("Z coordinate: ");
                        display_charge_properties(find_charge(find_charge_bycoord(x, y, z)));
                        input('Press ENTER to continue!');
                    end
                end

            case 5
                plot_ch();
                info_text;
                break
                
            case 7
                clear_space();

            otherwise
                error("Invalid Choice.");
        end
    end
end


% Function that creates the Calculate Values with Coulomb's Law menu
function calculator()
    global charge_space; global charge_space_permittivity; choice = 0;
    while choice ~= 6
        clc;
        fprintf("\t\t\t\t\t\t\t\t\tCALCULATE VALUES WITH COULOMB's LAW");
        if (charge_space(1).x_coord == 'N')
            n_charges = 0;
        else
            n_charges = length(charge_space);
        end
        fprintf("\n")
        disp(['The Charge space currently has ', num2str(n_charges), ' charges. The permittivity of the medium is ', num2str(charge_space_permittivity), '.']);
        disp("Enter 1 to change permittivity of Charge Space Medium: ");
        disp("Enter 2 to calculate superpositioned Net Force on any charge already present in Charge Space: ");
        disp("Enter 3 to calculate force for two custom charges, based on parameters entered for two charges: ");
        disp("Enter 4 to calculate force between any two charges already present in charge space: ");
        disp("Enter 5 to calculate superpositioned Electric field at any point in the Charge Space: ");
        disp("Enter 6 to exit menu: ");
        choice = input("==> ");
        clc;
        switch choice
            case 1
                ch = 0;
                while ch ~= 1 && ch ~= 2
                    disp("Enter 1 to manually enter Permittivity:");
                    disp("Enter 2 to multiply current permittivity by a relative permittivity: ");
                    ch = input("===>");
                    if ch == 1
                        fprintf('\n\n');
                        charge_space_permittivity = input("Enter new Permittivity of Charge Space: ");
                    elseif ch == 2
                        fprintf('\n\n');
                        charge_space_permittivity = charge_space_permittivity * input("Enter relative Permittivity of Charge Space: ");
                    end
                end   
                
            case 2
                ch = 0;
                while ch ~= 1 && ch ~= 2
                    if charge_space(1).x_coord == 'N'
                        disp('Charge Space is empty! Returning zero Force');
                        fx = 0;
                        fy = 0;
                        fz = 0;
                        fprintf('\nThe net force is:\n\t%e N along X-axis\n\t%e N along Y-axis\n\t%e N along Z-axis\n', fx, fy, fz);
                        input('Press ENTER to continue!');
                        break;
                    end
                    disp("Enter 1 to select Charge by it's number in Charge Space: ");
                    disp("Enter 2 to select Charge by  it's coordinate in Charge Space: ");
                    ch = input("===>");
                    fprintf('\n');
                    if ch == 1
                        rem_num = input("Enter the number of the charge on which you want to find net force: ");
                        [fx, fy, fz] = net_force_on(rem_num);
                        loco = ['charge ', num2str(rem_num), ' of Charge Space'];
                    elseif ch == 2
                        disp("Enter the x, y, z coordinates of the charge on which you want to find net force:"); 
                        x = input("X coordinate: ");
                        y = input("Y coordinate: ");
                        z = input("Z coordinate: ");
                        [fx, fy, fz] = net_force_on(find_charge_bycoord(x, y, z));
                        loco = ['charge at (', num2str(x), ', ', num2str(y), ', ', num2str(y), ')'];
                    end
                    fprintf('\nThe net force on %s is:\n\t%e N along X-axis\n\t%e N along Y-axis\n\t%e N along Z-axis\n', loco, fx, fy, fz);
                    input('Press ENTER to continue!')
                end

            case 3
                disp('Enter the x, y, z coordinates and the magnitude of the first charge: ');
                x = input("X coordinate: ");
                y = input("Y coordinate: ");
                z = input("Z coordinate: ");
                mag = input("Magnitude: ");
                loco1 = ['(', num2str(x), ', ', num2str(y), ', ', num2str(y), ')'];
                cA = charge_create(x, y, z, mag, 0);
                fprintf('\n');
                disp('Enter the x, y, z coordinates and the magnitude of the second charge: ');
                x = input("X coordinate: ");
                y = input("Y coordinate: ");
                z = input("Z coordinate: ");
                mag = input("Magnitude: ");
                loco2 = ['(', num2str(x), ', ', num2str(y), ', ', num2str(y), ')'];
                cB = charge_create(x, y, z, mag, 0);
                fprintf('\n');                   
                [fx, fy, fz] = force_on_two(cA, cB);
                fprintf('\nThe force on charge at %s by charge at %s is:\n\t%e N along X-axis\n\t%e N along Y-axis\n\t%e N along Z-axis\n', loco1, loco2, fx, fy, fz);
                input('Press ENTER to continue!')

            case 4
                if charge_space(1).x_coord == 'N'
                        disp('Charge Space is empty! Returning zero Force');
                        fx = 0;
                        fy = 0;
                        fz = 0;
                        fprintf('\nThe net force is:\n\t%e N along X-axis\n\t%e N along Y-axis\n\t%e N along Z-axis\n', fx, fy, fz);
                        input('Press Enter to continue.');
                else 
                    ch = 0;
                    while ch ~= 1 && ch ~= 2
                        disp("Choosing Charge on which Force will be measured:");
                        disp("Enter 1 to select it by number of Charge out of Charge Space: ");
                        disp("Enter 2 to select it by  it's coordinate in Charge Space: ");
                        ch = input("===>");
                        if ch == 1
                            rem_num = input("Enter the number of the charge on which you want to find force: ");
                            cA = rem_num;
                        elseif ch == 2
                            disp("Enter the x, y, z coordinates of the charge on which you want to find force:"); 
                            x = input("X coordinate: ");
                            y = input("Y coordinate: ");
                            z = input("Z coordinate: ");
                            cA = find_charge_bycoord(x, y, z);
                        end
                    end
                    ch = 0;
                    while ch ~= 1 && ch ~= 2
                        disp("Choosing Charge because of which Force will be experienced on first charge:");
                        disp("Enter 1 to select it by number of Charge out of Charge Space: ");
                        disp("Enter 2 to select it by  it's coordinate in Charge Space: ");
                        ch = input("===>");
                        if ch == 1
                            rem_num = input("Enter the number of the charge because of which Force is experienced: ");
                            cB = rem_num;
                        elseif ch == 2
                            disp("Enter the x, y, z coordinates of the charge because of which Force is experienced:"); 
                            x = input("X coordinate: ");
                            y = input("Y coordinate: ");
                            z = input("Z coordinate: ");
                            cB = find_charge_bycoord(x, y, z);
                        end
                    end
                    [fx, fy, fz] = force_on_two(cA, cB);
                    fprintf('\nThe force on first charge by second charge is:\n\t%e N along X-axis\n\t%e N along Y-axis\n\t%e N along Z-axis\n', fx, fy, fz);
                    input('Press ENTER to continue!')
                end
                
                
            case 5
                disp("Enter the x, y, z coordinates of the point where you want to find net Electric Field Strength:"); 
                x = input("X coordinate: ");
                y = input("Y coordinate: ");
                z = input("Z coordinate: ");
                [fx, fy, fz] = net_field_on(x, y, z);
                loco = ['(', num2str(x), ', ', num2str(y), ', ', num2str(y), ')'];
                fprintf('\nThe Electric Field at %s is:\n\t%e N/C along X-axis\n\t%e N/C along Y-axis\n\t%e N/C along Z-axis\n', loco, fx, fy, fz);
                input('Press ENTER to continue!')

            case 6
                info_text();
                break;

            otherwise
            error("Invalid Choice.");
        end
    end
end


function custom_grapher()
    clc;
    disp("This function tests Force vs Distance, Force vs Charge product and Force vs permittivities ");
    disp("While keeping every other parameter constant in each case to observe effects of the tested parameter.");
    fprintf("\n");
    low = input("Enter varied parameter's lower limit: ");
    high = input("Enter varied parameter's upper limit: ");
    step = input("Enter number of instances of the varied parameter: ");
    
    varied = linspace(low, high, step);
    charges = randi(high) * ones(1, step);
    distances = randi(high) * ones(1, step);
    permittivities = randi(high) * ones(1, step);
    
    % Scalar Distances are Varied
    forces = charges ./ (4 * pi * permittivities .* varied.^2);
    subplot(3, 1, 1);
    plot(varied, forces);
    xlabel("Scalar Distance between charge pairs");
    ylabel("Force acting on the Charge pair");
    title("Force vs Distance Curve");
    
    % Charge products are Varied
    forces = varied ./ (4 * pi * permittivities .* distances.^2);
    subplot(3, 1, 2);
    plot(varied, forces);
    xlabel("Product of Magnitude of Charges");
    ylabel("Force acting on the Charge pair");
    title("Force vs Charges' product Curve");
    
    % Permittivities are Varied
    forces = charges ./ (4 * pi * varied .* distances.^2);
    subplot(3, 1, 3);
    plot(varied, forces);
    xlabel("Permittivity of Medium");
    ylabel("Force acting on the Charge pair");
    title("Force vs Permittivity Curve");
    clc();
    info_text();
end


% Plots all charges in 3 Dimensional Charge Space as a Scatter Plot
function plot_ch()
    global charge_space;
    if charge_space(1).x_coord == 'N'
        disp("Empty Charge Space!")
        return
    end
    num = length(charge_space);
    x = zeros(1, num);
    y = zeros(1, num);
    z = zeros(1, num);
    sizes = zeros(1, num);
    colors = repmat([0, 0, 0], num, 1);
    for c = 1 : num
        charge = charge_space(c);
        x(c) = charge.x_coord;
        y(c) = charge.y_coord;
        z(c) = charge.z_coord;
        sizes(c) = charge.mag;
        if charge.col == 'r'
            colors(c, :) = [1, 0, 0];
        elseif charge.col == 'b'
            colors(c, :) = [0, 0, 1];
        end
    end
    sizes = (log(sizes) + 10^(-1)) * 25;
    sizes = abs(sizes);
    scatter3(x, y, z, sizes, colors, 'filled')
    title("Charge Space");
    xlabel("X axis Coordinate");
    ylabel("Y axis Coordinate");
    zlabel("Z axis Coordinate");
end


function [force_x, force_y, force_z] = net_force_on(obj)
    global charge_space; global charge_space_permittivity;
    flag = 0;
    
    % Finding the charge
    if isnumeric(obj)
        this_charge = find_charge(obj);
        if this_charge.x_coord == 'N'
            disp('Charge Space is empty! Returning zero Force');
            force_x = 0;
            force_y = 0;
            force_z = 0;
            return
        end
    else
        for charge = 1 : length(charge_space)
            if isequal(obj, charge_space(charge))
                flag = 1;
                this_charge = charge_space(charge);
                break
            end
        end
        if flag == 0
           disp('Charge Not Found!');
           return
        end
    end
    
    % If charge got, calculations to find the force and it's components
    vect_forces = zeros(3, length(charge_space) - 1);
    for c = 1 : length(charge_space)
        charge = charge_space(c);
        if isequal(charge, this_charge) || this_charge.x_coord == 'N'
            continue
        end
        charge_mag = charge.mag;
        charge_dist = dist_get(charge.x_coord, this_charge.x_coord, charge.y_coord, this_charge.y_coord, charge.z_coord, this_charge.z_coord);
        charge_dist = charge_dist ^ 3;
        vect_forces(:, c) = [
            this_charge.x_coord - charge.x_coord;
            this_charge.y_coord - charge.y_coord;
            this_charge.z_coord - charge.z_coord
        ];
        vect_forces(:, c) =  charge_mag / charge_dist * vect_forces(:, c);
    end
    comp_x = sum(vect_forces(1, :));
    comp_y = sum(vect_forces(2, :));
    comp_z = sum(vect_forces(3, :));
    force = this_charge.mag / (4 * pi * charge_space_permittivity);
    force_x = comp_x * force;
    force_y = comp_y * force;
    force_z = comp_z * force;
end


function [force_x, force_y, force_z] = force_on_two(obj, obj2)
    global charge_space_permittivity;
    
    % Finding the charges
    if isnumeric(obj)
        this_charge = find_charge(obj);
    else
        this_charge = obj;
    end
    
    if isnumeric(obj2)
        this_charge2 = find_charge(obj2);
    else
        this_charge2 = obj2;
    end
    
    % If charges found, calculations to find the force and it's components
    charge_dist = dist_get(this_charge.x_coord, this_charge2.x_coord, this_charge.y_coord, this_charge2.y_coord, this_charge.z_coord, this_charge2.z_coord);
    vect_dists = [
        this_charge.x_coord - this_charge2.x_coord;
        this_charge.y_coord - this_charge2.y_coord;
        this_charge.z_coord - this_charge2.z_coord
    ];
    charge_dist = charge_dist.^3;
    magnitude = (this_charge.mag * this_charge2.mag) / charge_dist;
    comp_x = sum(vect_dists(1));
    comp_y = sum(vect_dists(2));
    comp_z = sum(vect_dists(3));
    force = (magnitude) / (4 * pi * charge_space_permittivity);
    force_x = comp_x * force;
    force_y = comp_y * force;
    force_z = comp_z * force;
end


function [field_x, field_y, field_z] = net_field_on(x, y, z)
    global charge_space; global charge_space_permittivity;
   
    % Calculations to find the Electric Field's components
    vect_fields = zeros(3, length(charge_space) - 1);
    for c = 1 : length(charge_space)
        charge = charge_space(c);
        if charge.x_coord == 'N'
            continue
        end
        charge_mag = charge.mag;
        charge_dist = dist_get(charge.x_coord, x, charge.y_coord, y, charge.z_coord, z);
        charge_dist = charge_dist ^ 3;
        vect_fields(:, c) = [
            x - charge.x_coord;
            y - charge.y_coord;
            z - charge.z_coord
        ];
        vect_fields(:, c) =  charge_mag / charge_dist * vect_fields(:, c);
    end
    comp_x = sum(vect_fields(1, :));
    comp_y = sum(vect_fields(2, :));
    comp_z = sum(vect_fields(3, :));
    force = 1 / (4 * pi * charge_space_permittivity);
    field_x = comp_x * force;
    field_y = comp_y * force;
    field_z = comp_z * force;
end


function charge = find_charge(charge_number)
    global charge_space;
    if charge_number < 1 || charge_number > length(charge_space)
       disp('Enter a valid Charge number from Charge Space!');
    end
    charge = charge_space(charge_number);
end


function charge_num = find_charge_bycoord(x, y, z)
    global charge_space
    count = 1;
    for charge = charge_space
        if charge.x_coord == x && charge.y_coord == y && charge.z_coord == z
            charge_num = count;
            count = -1;
            break
        end
        count = count + 1;
    end
    if count ~= -1
        disp("No charge at that exact Coordinate!")
        charge_num = 0;
    end
end


function charge_dest(obj, show)
    global charge_space
    if nargin == 1
        show = 1;
    end
    flag = 0;
    if isnumeric(obj)
        if obj < 1 || obj > length(charge_space)
            disp('Enter a valid Charge number from Charge Space!');
            input("");
            return
        elseif charge_space(1).x_coord == 'N'
            disp("Charge space is empty. No charges to remove!");
            input("");
            return
        end
        charge_space(obj) = [];
    else
        for charge = 1 : length(charge_space)
            if isequal(obj, charge_space(charge))
                flag = 1;
                if length(charge_space) == 1
                    charge_space(charge).x_coord = 'N';
                end
                charge_space(charge) = [];
                break
            end
        end
        if flag == 0
           disp('Charge Not Found!');
        end
    end
    if show == 1
        plot_ch()
    end
end


function charge_out = charge_create(x_coordinate, y_coordinate, z_coordinate, magnitude, show)
    global charge_space;
    if nargin == 4
        show = 1;
    end
    if magnitude == 0
        disp('Zero Magnitude Charge created! It has no effect on Charge Space and will be ignored.');
        charge_out.x_coord = x_coordinate;
        charge_out.y_coord = y_coordinate;  
        charge_out.z_coord = z_coordinate;
        charge_out.mag = magnitude;
        return
    end
    for charge = 1 : length(charge_space)
        if x_coordinate == charge_space(charge).x_coord && y_coordinate == charge_space(charge).y_coord && z_coordinate == charge_space(charge).z_coord
            charge_space(charge).mag = magnitude;
            if magnitude > 0
                charge_space(charge).col = 'r';
            elseif magnitude == 0
                charge_space(charge).col = 'g';
            else
                charge_space(charge).col = 'b';
            end
            charge_out = charge_space(charge);
            disp('Charge Already exists here in charge space. Magnitude updated.');
            if show == 1
                plot_ch()
            end
            return;
        end
    end
    charge_out.x_coord = x_coordinate;
    charge_out.y_coord = y_coordinate;  
    charge_out.z_coord = z_coordinate;
    charge_out.mag = magnitude;
    if magnitude > 0
        charge_out.col = 'r';
    else
        charge_out.col = 'b';
    end
    if charge_space(end).x_coord == 'N'
        charge_space(end) = charge_out;
    else
        charge_space(end + 1) = charge_out;
    end
    if show == 1
        plot_ch()
    end
end


function clear_space()
    global charge_space
    charge_space = struct('x_coord', 'N', 'y_coord', [], 'z_coord', [], 'mag', [], 'col', 'g');
end


function create_rnd_charges(number)
    for i = 1: number
        charge_create(randi([-10000, 10000]), randi([-10000, 10000]), randi([-10000, 10000]), randi([-1000, 1000]), 0);
    end
end


function display_charge_properties(charge)
    fprintf('\n');
    if charge.x_coord == 'N'
        disp("Charge Space is empty! Make sure you place charges! ");
        return
    end
    disp(['Location of this charge in 3d space is (', ...
        num2str(charge.x_coord), ', ', ...
        num2str(charge.y_coord), ', ', ...
        num2str(charge.z_coord), ').']);
    disp(['Magnitude of the charge is ', num2str(charge.mag), '.']);
end


function info_text()
    global charge_space; global charge_space_permittivity;
    if (charge_space(1).x_coord == 'N')
        n_charges = 0;
    else
        n_charges = length(charge_space);
    end
    disp(['The Charge space currently has ', num2str(n_charges), ' charges. The permittivity of the medium is ', num2str(charge_space_permittivity), '.']);
    fprintf("Call menu() function to add and remove charges in Charge Space or plot Charge Space in three dimensions.");
    fprintf("\nCall calc() function to calculate Forces, Field after adding charges in Charge Space.");
    fprintf("\nCall plot_graphs() function to plot Force vs Distance or Charge curves with custom Charge pairs.\n"); 
end


function distance = dist_get(x1, x2, y1, y2, z1, z2)
    distance = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
end