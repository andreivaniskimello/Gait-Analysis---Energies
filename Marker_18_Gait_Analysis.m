%18 Points Mechanical Energy Routine
%by André Ivaniski Mello (andreivaniskimello@gmail.com)


%INPUT: 
% Input File and Output File Paths
% Subject Name and Condition

% The INPUT FILE is a .txt file with kinematic data of 18 markers points positions through time; 

% INSERT ALSO:

% Touch down (TD) and Take Off (TO) in FRAMES of STRIDES to be analyzed

% Treadmill Speed, Treadmill Inclination, Body Mass
% Sample Frequency
% Butterworth Filter Parameters desired


%OUTPUT:

%Excel file with

% Spatiotemporal parameters of locomotion per Stride
% Mechanical Work per Stride
% Mechanical Power per Stride


%X Axis: ML
%Y Axis: AP
%Z Axis: Vertical

clear all
close all
clc

% THIS ROUTINE HAS 09 SECTIONS

% SECTION 01: INFORMATION FOR ROUTINE 
% SECTION 02: IMPORT DATA AND INFORMATION PREPARATION 
% SECTION 03: LINEAR SPATIOTEMPORAL ANALYSIS
% SECTION 04: CENTER OF MASS POSITIONS
% SECTION 05: ANGULAR ANALYSIS
% SECTION 06: ENERGIES
% SECTION 07: TOTAL MECHANICAL ENERGY
% SECTION 08: EXPORT DATA
% SECTION 09: GRAPHICS

%**************************************************************************
%% SECTION 01
% INFORMATION FOR ROUTINE %

%Insert 'SubjectName_Speed' here, in order to export the Output Files
%properly

File_Name =['Herculano_80'];  

%Insert Input File Path here
Kinematic_file = ('C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto TPAMA - Leonardo Bloedow\Dados\Herculano\HERCULANO 80_Edited.txt');

%Insert Output File Path here
Output_File_Path = ['C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto TPAMA - Leonardo Bloedow\Dados\Data Out'];     

%Subject Details
Body_Mass = 80;               %in kg
Treadmill_Speed = 5.3056;     %insert treadmill speed in m/s


%STRIDE 01       
St1_TD1 = 2358;
St1_TO = 2387;
St1_TD2 = 2460;

%STRIDE 2
St2_TD1 = 2460;
St2_TO = 2492;
St2_TD2 = 2564;

%STRIDE 3
St3_TD1 = 2564;
St3_TO = 2599;
St3_TD2 = 2670;

%% SECTION 02
% IMPORT DATA AND INFORMATION PREPARATION %

g = 9.81;
Sample_Frequency = 200;                                  %in Hertz
dt = 1/Sample_Frequency;
Treadmill_Inclination = 0.25;                            %insert inclination as decimal. Ex: if inclination is "25%", insert as "0.25"
Vt = Treadmill_Speed;                                    %Vt: velocity of treadmill
Treadmill_Angle = atand(Treadmill_Inclination);          %Inclination of treadmill in degrees
Vt_Y = Vt * cosd(Treadmill_Angle);
Vt_Z = Vt * sind(Treadmill_Angle);

% Load kinematic data

delimiterIn = '\t';
Kinematic = importdata(Kinematic_file,delimiterIn);

%Stride Analyzed
T0 = Kinematic(1,1);

Frames = size(Kinematic,1);                              % Calculate the total frames number
Time_Total = Frames/Sample_Frequency;                    %Time = total frames / fsample
Time = (0:dt:Time_Total);
Time_Row = Time.';


St1_TD1 = St1_TD1 - T0;
St1_TO = St1_TO - T0;
St1_TD2 = St1_TD2 - T0;

St2_TD1 = St2_TD1 - T0;
St2_TO = St2_TO - T0;
St2_TD2 = St2_TD2 - T0;

St3_TD1 = St3_TD1 - T0;
St3_TO = St3_TO - T0;
St3_TD2 = St3_TD2 - T0;

Kinematic = Kinematic.*0.001;                           %milimeters (mm) to meters (m).

fcut = 20;
order = 2;
Kinematic = matfiltfilt(dt, fcut, order, Kinematic);

R_Wrist = Kinematic(:,3:5);                 %Right_Wrist
R_Elbow = Kinematic(:,6:8);                 %Right_Elbow
R_Shoulder = Kinematic(:,9:11);             %Right_Shoulder

L_Wrist = Kinematic(:,12:14);               %Left_Wrist
L_Elbow = Kinematic(:,15:17);               %Left_Elbow
L_Shoulder = Kinematic(:,18:20);            %Left_Shoulder

R_Throc = Kinematic(:,21:23);               %Right_Throcanter
R_Knee = Kinematic(:,24:26);                %Right_Knee
R_Ankle = Kinematic(:,27:29);               %Right_Ankle
R_Heel = Kinematic(:,30:32);                %Right_Heel
R_Toe = Kinematic(:,33:35);                 %Right_Toe

L_Throc = Kinematic(:,36:38);               %Left_Throcanter
L_Knee = Kinematic(:,39:41);                %Left_Knee
L_Ankle = Kinematic(:,42:44);               %Left_Ankle
L_Heel = Kinematic(:,45:47);                %Left_Heel
L_Toe = Kinematic(:,48:50);                 %Left_Toe

R_Head = Kinematic(:,51:53);                %Right_Head
L_Head = Kinematic(:,54:56);                %Left_Head


Medium_Throc = (R_Throc + L_Throc)/2;
Medium_Head = (R_Head + L_Head)/2;

%% SECTION 03
% LINEAR SPATIOTEMPORAL ANALYSIS %

       %STRIDE 01
St1_Duration = (St1_TD2 - St1_TD1)*dt;
St1_Frequency = 1/St1_Duration;
St1_Lenght = Vt/St1_Frequency;
St1_Speed = St1_Frequency * St1_Lenght;
St1_Contact_Phase_Duration = (St1_TO - St1_TD1)*dt;
St1_Duty_Factor = (St1_Contact_Phase_Duration/St1_Duration)*100;

St1_Time_Vector = [0:dt:((St1_TD2-St1_TD1)/Sample_Frequency)];
St1_Time_Stride_Percentage = St1_Time_Vector/(max(St1_Time_Vector))*100; 


       %STRIDE 02
St2_Duration = (St2_TD2 - St2_TD1)*dt;
St2_Frequency = 1/St2_Duration;
St2_Lenght = Vt/St2_Frequency;
St2_Speed = St2_Frequency * St2_Lenght;
St2_Contact_Phase_Duration = (St2_TO - St2_TD1)*dt;
St2_Duty_Factor = (St2_Contact_Phase_Duration/St2_Duration)*100;

St2_Time_Vector = [0:dt:((St2_TD2-St2_TD1)/Sample_Frequency)];
St2_Time_Stride_Percentage = St2_Time_Vector/(max(St2_Time_Vector))*100;


       %STRIDE 03
St3_Duration = (St3_TD2 - St3_TD1)*dt;
St3_Frequency = 1/St3_Duration;
St3_Lenght = Vt/St3_Frequency;
St3_Speed = St3_Frequency * St3_Lenght;
St3_Contact_Phase_Duration = (St3_TO - St3_TD1)*dt;
St3_Duty_Factor = (St3_Contact_Phase_Duration/St3_Duration)*100;

St3_Time_Vector = [0:dt:((St3_TD2-St3_TD1)/Sample_Frequency)];
St3_Time_Stride_Percentage = St3_Time_Vector/(max(St3_Time_Vector))*100;


%% SECTION 04
% CENTER OF MASS POSITIONS (per segment and all body) %

L_Foot_CoM = ((L_Toe - L_Ankle) *0.50 + L_Ankle);           % Foot Center of Mass
L_Shank_CoM = ((L_Ankle - L_Knee)*0.433 + L_Knee);          % Shank Center of Mass
L_Thigh_CoM = ((L_Knee - L_Throc)*0.433 + L_Throc);         % Thigh Center of Mass
L_Trunk_CoM = ((L_Throc - L_Head)*0.66 + L_Head);           % Trunk Center of Mass
L_Arm_CoM = ((L_Shoulder - L_Elbow)*0.436 + L_Elbow);       % Arm Center of Mass
L_Forearm_CoM = ((L_Elbow - L_Wrist)*0.682 + L_Wrist);      % Forearm and Hand Center of Mass

R_Foot_CoM = ((R_Toe - R_Ankle) *0.50 + R_Ankle);           % Foot Center of Mass
R_Shank_CoM = ((R_Ankle - R_Knee)*0.433 + R_Knee);          % Shank Center of Mass
R_Thigh_CoM = ((R_Knee - R_Throc)*0.433 + R_Throc);         % Thigh Center of Mass
R_Trunk_CoM = ((R_Throc - R_Head)*0.66 + R_Head);           % Trunk Center of Mass
R_Arm_CoM = ((R_Shoulder - R_Elbow)*0.436 + R_Elbow);       % Arm Center of Mass
R_Forearm_CoM = ((R_Elbow - R_Wrist)*0.682 + R_Wrist);      % Forearm and Hand Center of Mass

Medium_Trunk_CoM = (L_Trunk_CoM + R_Trunk_CoM)/2;

Foot_Mass= 0.0145*Body_Mass;
Shank_Mass = 0.0465*Body_Mass;
Thigh_Mass = 0.100*Body_Mass;
Trunk_Mass = 0.289*Body_Mass;                               %Half trunk
Trunk_Mass = Trunk_Mass*2;
Arm_Mass = 0.028*Body_Mass;
Forearm_Hand_Mass = 0.022*Body_Mass;


% BODY CENTER OF MASS
CoM(:,1:3) = ((R_Foot_CoM).*Foot_Mass)+((R_Shank_CoM).*Shank_Mass)+((R_Thigh_CoM).*Thigh_Mass)+((R_Trunk_CoM).*Trunk_Mass)+((R_Arm_CoM).*Arm_Mass)+((R_Forearm_CoM).*Forearm_Hand_Mass)+ ((L_Foot_CoM).*Foot_Mass)+((L_Shank_CoM).*Shank_Mass)+((L_Thigh_CoM).*Thigh_Mass)+((L_Trunk_CoM.*Trunk_Mass)+((L_Arm_CoM).*Arm_Mass)+((L_Forearm_CoM).*Forearm_Hand_Mass));

[nl, nc]=size(CoM);

CoM(:,:) = CoM(:,:)/100;
CoM_X = CoM(:,1);
CoM_Y = CoM(:,2);
CoM_Z = CoM(:,3);

fcut_d = 8;
CoM_Y = matfiltfilt(dt, fcut_d, order, CoM_Y);
CoM_X = matfiltfilt(dt, fcut_d, order, CoM_X);
CoM_Z = matfiltfilt(dt, fcut_d, order, CoM_Z);


%CoM  VERTICAL POSITION CORRECTION BY THE TREADMILL INCLINATION

%Vertical velocity (Vv)
Vv = Vt*(sind(Treadmill_Angle));

%Vertical displacement (Sv)
Sv = Time*Vv;
Sv = Sv.';
Sv(1) = [];

% CoM Vertical Axis Corrected by the Inclination
CoM_Slope = CoM(:,3)+ Sv;
  
%% SECTION 05
% ANGULAR ANALYSIS %

% Segments Vector Determination

R_Foot = R_Toe - R_Heel;

R_Foot(:,1) = (R_Foot(:,1)-(min(R_Foot(:,1))))+(abs(min(R_Foot(:,1))));
R_Foot(:,2) = (R_Foot(:,2)-(min(R_Foot(:,2))))+(abs(min(R_Foot(:,2))));
R_Foot(:,3) = (R_Foot(:,3)-(min(R_Foot(:,3))))+(abs(min(R_Foot(:,3))));

R_Shank = R_Ankle - R_Knee;

R_Shank(:,1) = (R_Shank(:,1)-(min(R_Shank(:,1))))+(abs(min(R_Shank(:,1))));
R_Shank(:,2) = (R_Shank(:,2)-(min(R_Shank(:,2))))+(abs(min(R_Shank(:,2))));
R_Shank(:,3) = (R_Shank(:,3)-(min(R_Shank(:,3))))+(abs(min(R_Shank(:,3))));

R_Thigh = R_Knee - R_Throc;

R_Thigh(:,1) = (R_Thigh(:,1)-(min(R_Thigh(:,1))))+(abs(min(R_Thigh(:,1))));
R_Thigh(:,2) = (R_Thigh(:,2)-(min(R_Thigh(:,2))))+(abs(min(R_Thigh(:,2))));
R_Thigh(:,3) = (R_Thigh(:,3)-(min(R_Thigh(:,3))))+(abs(min(R_Thigh(:,3))));

R_Arm = R_Shoulder - R_Elbow;

R_Arm(:,1) = (R_Arm(:,1)-(min(R_Arm(:,1))))+(abs(min(R_Arm(:,1))));
R_Arm(:,2) = (R_Arm(:,2)-(min(R_Arm(:,2))))+(abs(min(R_Arm(:,2))));
R_Arm(:,3) = (R_Arm(:,3)-(min(R_Arm(:,3))))+(abs(min(R_Arm(:,3))));




R_Forearm = R_Elbow - R_Wrist;

R_Forearm(:,1) = (R_Forearm(:,1)-(min(R_Forearm(:,1))))+(abs(min(R_Forearm(:,1))));
R_Forearm(:,2) = (R_Forearm(:,2)-(min(R_Forearm(:,2))))+(abs(min(R_Forearm(:,2))));
R_Forearm(:,3) = (R_Forearm(:,3)-(min(R_Forearm(:,3))))+(abs(min(R_Forearm(:,3))));

L_Foot = L_Toe - L_Heel;

L_Foot(:,1) = (L_Foot(:,1)-(min(L_Foot(:,1))))+(abs(min(L_Foot(:,1))));
L_Foot(:,2) = (L_Foot(:,2)-(min(L_Foot(:,2))))+(abs(min(L_Foot(:,2))));
L_Foot(:,3) = (L_Foot(:,3)-(min(L_Foot(:,3))))+(abs(min(L_Foot(:,3))));

L_Shank = L_Ankle - L_Knee;

L_Shank(:,1) = (L_Shank(:,1)-(min(L_Shank(:,1))))+(abs(min(L_Shank(:,1))));
L_Shank(:,2) = (L_Shank(:,2)-(min(L_Shank(:,2))))+(abs(min(L_Shank(:,2))));
L_Shank(:,3) = (L_Shank(:,3)-(min(L_Shank(:,3))))+(abs(min(L_Shank(:,3))));

L_Thigh = L_Knee - L_Throc;

L_Thigh(:,1) = (L_Thigh(:,1)-(min(L_Thigh(:,1))))+(abs(min(L_Thigh(:,1))));
L_Thigh(:,2) = (L_Thigh(:,2)-(min(L_Thigh(:,2))))+(abs(min(L_Thigh(:,2))));
L_Thigh(:,3) = (L_Thigh(:,3)-(min(L_Thigh(:,3))))+(abs(min(L_Thigh(:,3))));

L_Arm = L_Shoulder - L_Elbow;

L_Arm(:,1) = (L_Arm(:,1)-(min(L_Arm(:,1))))+(abs(min(L_Arm(:,1))));
L_Arm(:,2) = (L_Arm(:,2)-(min(L_Arm(:,2))))+(abs(min(L_Arm(:,2))));
L_Arm(:,3) = (L_Arm(:,3)-(min(L_Arm(:,3))))+(abs(min(L_Arm(:,3))));

L_Forearm = L_Elbow - L_Wrist;

L_Forearm(:,1) = (L_Forearm(:,1)-(min(L_Forearm(:,1))))+(abs(min(L_Forearm(:,1))));
L_Forearm(:,2) = (L_Forearm(:,2)-(min(L_Forearm(:,2))))+(abs(min(L_Forearm(:,2))));
L_Forearm(:,3) = (L_Forearm(:,3)-(min(L_Forearm(:,3))))+(abs(min(L_Forearm(:,3))));

Medium_Trunk = (Medium_Throc - Medium_Head);


     %%% SEGMENTS ANGLES IN SAGITAL PLANE (YZ) %%%

R_Foot_Angle(:,1) = atan(R_Foot(:,2)./R_Foot(:,3));
R_Shank_Angle(:,1) = atan(R_Shank(:,2)./R_Shank(:,3));
R_Thigh_Angle(:,1) = atan(R_Thigh(:,2)./R_Thigh(:,3));

R_Arm_Angle(:,1) = atan(R_Arm(:,2)./R_Arm(:,3));
R_Forearm_Angle(:,1) = atan(R_Forearm(:,2)./R_Forearm(:,3));

L_Foot_Angle(:,1) = atan(L_Foot(:,2)./L_Foot(:,3));
L_Shank_Angle(:,1) = atan(L_Shank(:,2)./L_Shank(:,3));
L_Thigh_Angle(:,1) = atan(L_Thigh(:,2)./L_Thigh(:,3));

L_Arm_Angle(:,1) = atan(L_Arm(:,2)./L_Arm(:,3));
L_Forearm_Angle(:,1) = atan(L_Forearm(:,2)./L_Forearm(:,3));

Medium_Trunk_Angle(:,1) = atan(Medium_Trunk(:,2)./Medium_Trunk(:,3));


       %%% SEGMENTS ANGLES IN FRONTAL PLANE (XZ) %%%

R_Foot_Angle(:,2) = atan(R_Foot(:,1)./R_Foot(:,3));
R_Shank_Angle(:,2) = atan(R_Shank(:,1)./R_Shank(:,3));
R_Thigh_Angle(:,2) = atan(R_Thigh(:,1)./R_Thigh(:,3));

R_Arm_Angle(:,2) = atan(R_Arm(:,1)./R_Arm(:,3));
R_Forearm_Angle(:,2) = atan(R_Forearm(:,1)./R_Forearm(:,3));

L_Foot_Angle(:,2) = atan(L_Foot(:,1)./L_Foot(:,3));
L_Shank_Angle(:,2) = atan(L_Shank(:,1)./L_Shank(:,3));
L_Thigh_Angle(:,2) = atan(L_Thigh(:,1)./L_Thigh(:,3));

L_Arm_Angle(:,2) = atan(L_Arm(:,1)./L_Arm(:,3));
L_Forearm_Angle(:,2) = atan(L_Forearm(:,1)./L_Forearm(:,3));

Medium_Trunk_Angle(:,2) = atan(Medium_Trunk(:,1)./Medium_Trunk(:,3));


    %%% SEGMENTS ANGLES IN TRANSVERSAL PLANE (XY) %%%

R_Foot_Angle(:,3) = atan(R_Foot(:,1)./R_Foot(:,2));
R_Shank_Angle(:,3) = atan(R_Shank(:,1)./R_Shank(:,2));
R_Thigh_Angle(:,3) = atan(R_Thigh(:,1)./R_Thigh(:,2));

R_Arm_Angle(:,3) = atan(R_Arm(:,1)./R_Arm(:,2));
R_Forearm_Angle(:,3) = atan(R_Forearm(:,1)./R_Forearm(:,2));

L_Foot_Angle(:,3) = atan(L_Foot(:,1)./L_Foot(:,2));
L_Shank_Angle(:,3) = atan(L_Shank(:,1)./L_Shank(:,2));
L_Thigh_Angle(:,3) = atan(L_Thigh(:,1)./L_Thigh(:,2));

L_Arm_Angle(:,3) = atan(L_Arm(:,1)./L_Arm(:,2));
L_Forearm_Angle(:,3) = atan(L_Forearm(:,1)./L_Forearm(:,2));

Medium_Trunk_Angle(:,3) = atan(Medium_Trunk(:,1)./Medium_Trunk(:,2));


% SEGMENTS ANGULAR SPEED (w) %

R_Foot_Angular_Speed = diff(R_Foot_Angle)./dt;
R_Shank_Angular_Speed = diff(R_Shank_Angle)./dt;
R_Thigh_Angular_Speed = diff(R_Thigh_Angle)./dt;
R_Arm_Angular_Speed = diff(R_Arm_Angle)./dt;
R_Forearm_Angular_Speed = diff(R_Forearm_Angle)./dt;

R_Foot_Angular_Speed = interpft(R_Foot_Angular_Speed,nl);
R_Shank_Angular_Speed = interpft(R_Shank_Angular_Speed,nl);
R_Thigh_Angular_Speed = interpft(R_Thigh_Angular_Speed,nl);
R_Arm_Angular_Speed = interpft(R_Arm_Angular_Speed,nl);
R_Forearm_Angular_Speed = interpft(R_Forearm_Angular_Speed,nl);

L_Foot_Angular_Speed = diff(L_Foot_Angle)./dt;
L_Shank_Angular_Speed = diff(L_Shank_Angle)./dt;
L_Thigh_Angular_Speed = diff(L_Thigh_Angle)./dt;
L_Arm_Angular_Speed = diff(L_Arm_Angle)./dt;
L_Forearm_Angular_Speed = diff(L_Forearm_Angle)./dt;
 
L_Foot_Angular_Speed = interpft(L_Foot_Angular_Speed,nl);
L_Shank_Angular_Speed = interpft(L_Shank_Angular_Speed,nl);
L_Thigh_Angular_Speed = interpft(L_Thigh_Angular_Speed,nl);
L_Arm_Angular_Speed = interpft(L_Arm_Angular_Speed,nl);
L_Forearm_Angular_Speed = interpft(L_Forearm_Angular_Speed,nl);

Medium_Trunk_Angular_Speed = diff(Medium_Trunk_Angle)./dt;
Medium_Trunk_Angular_Speed = interpft(Medium_Trunk_Angular_Speed,nl);        


%% SECTION 06
% ENERGIES %

        % POTENTIAL GRAVITATION ENERGY (Pg) (Pg = mgh)
        
Pg = Body_Mass * g * CoM_Slope;


        % KINETIC ENERGY OF BODY CENTER OF MASS (EK) (EK = m . v^2 * 1/2)
V_CoM(:,1:3) = diff(CoM)./dt;
V_CoM = interpft(V_CoM,nl);
V_CoM(:,2) = V_CoM(:,2) + Vt_Y;
V_CoM(:,3) = V_CoM(:,3) + Vt_Z;

EK_X = Body_Mass * (V_CoM(:,1).^2) * 0.5;
EK_Y = Body_Mass * (V_CoM(:,2).^2) * 0.5;
EK_Z = Body_Mass * (V_CoM(:,3).^2) * 0.5;
EK_Total = EK_X + EK_Y + EK_Z;

[nl nc] = size(EK_Total);



    %INTERNAL KINETIC ENERGY
                
        
         %Linear Velocity in ML (X Axis) direction of Each Segment Center of Mass
                    
    % Foot Right
V_R_Foot_CoM_Raw(:,1) = diff(R_Foot_CoM(:,1))./dt;                       % This first line calculate the Raw Linear Velocity of Segments Center of Mass
V_R_Foot_CoM_Interpft(:,1) = (interpft(V_R_Foot_CoM_Raw(:,1),nl));       % This second line interpolate the Raw Linear Velocity to the same amount of Lines of BOdy Center of Mass Velocity
V_R_Foot_CoM(:,1) = V_R_Foot_CoM_Interpft(:,1) - V_CoM(:,1);             % This third line subtracts the Body Center of Mass from the Segment Center of Mass

    % Shank Right
V_R_Shank_CoM_Raw(:,1) = diff(R_Shank_CoM(:,1))./dt; 
V_R_Shank_CoM_Interpft(:,1) = (interpft(V_R_Shank_CoM_Raw(:,1),nl));
V_R_Shank_CoM(:,1) = V_R_Shank_CoM_Interpft(:,1) - V_CoM(:,1);

    % Thigh Right
V_R_Thigh_CoM_Raw(:,1) = diff(R_Thigh_CoM(:,1))./dt; 
V_R_Thigh_CoM_Interpft(:,1) = (interpft(V_R_Thigh_CoM_Raw(:,1),nl));
V_R_Thigh_CoM(:,1) = V_R_Thigh_CoM_Interpft(:,1) - V_CoM(:,1);

    % Arm Right
V_R_Arm_CoM_Raw(:,1) = diff(R_Arm_CoM(:,1))./dt; 
V_R_Arm_CoM_Interpft(:,1) = (interpft(V_R_Arm_CoM_Raw(:,1),nl));
V_R_Arm_CoM(:,1) = V_R_Arm_CoM_Interpft(:,1) - V_CoM(:,1);

    % Forearm Right
V_R_Forearm_CoM_Raw(:,1) = diff(R_Forearm_CoM(:,1))./dt; 
V_R_Forearm_CoM_Interpft(:,1) = (interpft(V_R_Forearm_CoM_Raw(:,1),nl));
V_R_Forearm_CoM(:,1) = V_R_Forearm_CoM_Interpft(:,1) - V_CoM(:,1);


    % Foot Left
V_L_Foot_CoM_Raw(:,1) = diff(L_Foot_CoM(:,1))./dt;                      
V_L_Foot_CoM_Interpft(:,1) = (interpft(V_L_Foot_CoM_Raw(:,1),nl));       
V_L_Foot_CoM(:,1) = V_L_Foot_CoM_Interpft(:,1) - V_CoM(:,1);        

    % Shank Left
V_L_Shank_CoM_Raw(:,1) = diff(L_Shank_CoM(:,1))./dt; 
V_L_Shank_CoM_Interpft(:,1) = (interpft(V_L_Shank_CoM_Raw(:,1),nl));
V_L_Shank_CoM(:,1) = V_L_Shank_CoM_Interpft(:,1) - V_CoM(:,1);

    % Thigh Left
V_L_Thigh_CoM_Raw(:,1) = diff(L_Thigh_CoM(:,1))./dt; 
V_L_Thigh_CoM_Interpft(:,1) = (interpft(V_L_Thigh_CoM_Raw(:,1),nl));
V_L_Thigh_CoM(:,1) = V_L_Thigh_CoM_Interpft(:,1) - V_CoM(:,1);

    % Arm Left
V_L_Arm_CoM_Raw(:,1) = diff(L_Arm_CoM(:,1))./dt; 
V_L_Arm_CoM_Interpft(:,1) = (interpft(V_L_Arm_CoM_Raw(:,1),nl));
V_L_Arm_CoM(:,1) = V_L_Arm_CoM_Interpft(:,1) - V_CoM(:,1);

    % Forearm Left
V_L_Forearm_CoM_Raw(:,1) = diff(L_Forearm_CoM(:,1))./dt; 
V_L_Forearm_CoM_Interpft(:,1) = (interpft(V_L_Forearm_CoM_Raw(:,1),nl));
V_L_Forearm_CoM(:,1) = V_L_Forearm_CoM_Interpft(:,1) - V_CoM(:,1);

    % Medium Trunk
V_Trunk_Medium_CoM_Raw(:,1) = diff(Medium_Trunk_CoM(:,1))./dt; 
V_Trunk_Medium_CoM_Interpft(:,1) = (interpft(V_Trunk_Medium_CoM_Raw(:,1),nl));
V_Trunk_Medium_CoM(:,1) = V_Trunk_Medium_CoM_Interpft(:,1) - V_CoM(:,1);

    
        %Linear Velocity in AP (Y Axis) direction of Each Segment Center of Mass
    
    % Foot Right
V_R_Foot_CoM_Raw(:,2) = diff(R_Foot_CoM(:,2))./dt;  
V_R_Foot_CoM_Interpft(:,2) = (interpft(V_R_Foot_CoM_Raw(:,2),nl));
V_R_Foot_CoM(:,2) = V_R_Foot_CoM_Interpft(:,2) - V_CoM(:,2);

    % Shank Right
V_R_Shank_CoM_Raw(:,2) = diff(R_Shank_CoM(:,2))./dt; 
V_R_Shank_CoM_Interpft(:,2) = (interpft(V_R_Shank_CoM_Raw(:,2),nl));
V_R_Shank_CoM(:,2) = V_R_Shank_CoM_Interpft(:,2) - V_CoM(:,2);

    % Thigh Right
V_R_Thigh_CoM_Raw(:,2) = diff(R_Thigh_CoM(:,2))./dt; 
V_R_Thigh_CoM_Interpft(:,2) = (interpft(V_R_Thigh_CoM_Raw(:,2),nl));
V_R_Thigh_CoM(:,2) = V_R_Thigh_CoM_Interpft(:,2) - V_CoM(:,2);

    % Arm Right
V_R_Arm_CoM_Raw(:,2) = diff(R_Arm_CoM(:,2))./dt; 
V_R_Arm_CoM_Interpft(:,2) = (interpft(V_R_Arm_CoM_Raw(:,2),nl));
V_R_Arm_CoM(:,2) = V_R_Arm_CoM_Interpft(:,2) - V_CoM(:,2);

    % Forearm Right
V_R_Forearm_CoM_Raw(:,2) = diff(R_Forearm_CoM(:,2))./dt; 
V_R_Forearm_CoM_Interpft(:,2) = (interpft(V_R_Forearm_CoM_Raw(:,2),nl));
V_R_Forearm_CoM(:,2) = V_R_Forearm_CoM_Interpft(:,2) - V_CoM(:,2);


    % Foot Left
V_L_Foot_CoM_Raw(:,2) = diff(L_Foot_CoM(:,2))./dt;                      
V_L_Foot_CoM_Interpft(:,2) = (interpft(V_L_Foot_CoM_Raw(:,2),nl));       
V_L_Foot_CoM(:,2) = V_L_Foot_CoM_Interpft(:,2) - V_CoM(:,2);        

    % Shank Left
V_L_Shank_CoM_Raw(:,2) = diff(L_Shank_CoM(:,2))./dt; 
V_L_Shank_CoM_Interpft(:,2) = (interpft(V_L_Shank_CoM_Raw(:,2),nl));
V_L_Shank_CoM(:,2) = V_L_Shank_CoM_Interpft(:,2) - V_CoM(:,2);

    % Thigh Left
V_L_Thigh_CoM_Raw(:,2) = diff(L_Thigh_CoM(:,2))./dt; 
V_L_Thigh_CoM_Interpft(:,2) = (interpft(V_L_Thigh_CoM_Raw(:,2),nl));
V_L_Thigh_CoM(:,2) = V_L_Thigh_CoM_Interpft(:,2) - V_CoM(:,2);

    % Arm Left
V_L_Arm_CoM_Raw(:,2) = diff(L_Arm_CoM(:,2))./dt; 
V_L_Arm_CoM_Interpft(:,2) = (interpft(V_L_Arm_CoM_Raw(:,2),nl));
V_L_Arm_CoM(:,2) = V_L_Arm_CoM_Interpft(:,2) - V_CoM(:,2);

    % Forearm Left
V_L_Forearm_CoM_Raw(:,2) = diff(L_Forearm_CoM(:,2))./dt; 
V_L_Forearm_CoM_Interpft(:,2) = (interpft(V_L_Forearm_CoM_Raw(:,2),nl));
V_L_Forearm_CoM(:,2) = V_L_Forearm_CoM_Interpft(:,2) - V_CoM(:,2);

   % Medium Trunk
V_Trunk_Medium_CoM_Raw(:,2) = diff(Medium_Trunk_CoM(:,2))./dt; 
V_Trunk_Medium_CoM_Interpft(:,2) = (interpft(V_Trunk_Medium_CoM_Raw(:,2),nl));
V_Trunk_Medium_CoM(:,2) = V_Trunk_Medium_CoM_Interpft(:,2) - V_CoM(:,2);
  

        %Linear Velocity in V direction of Each Segment Center of Mass
    
     % Foot Right
V_R_Foot_CoM_Raw(:,3) = diff(R_Foot_CoM(:,3))./dt;  
V_R_Foot_CoM_Interpft(:,3) = (interpft(V_R_Foot_CoM_Raw(:,3),nl));
V_R_Foot_CoM(:,3) = V_R_Foot_CoM_Interpft(:,3) - V_CoM(:,3);

    % Shank Right
V_R_Shank_CoM_Raw(:,3) = diff(R_Shank_CoM(:,3))./dt; 
V_R_Shank_CoM_Interpft(:,3) = (interpft(V_R_Shank_CoM_Raw(:,3),nl));
V_R_Shank_CoM(:,3) = V_R_Shank_CoM_Interpft(:,3) - V_CoM(:,3);

    % Thigh Right
V_R_Thigh_CoM_Raw(:,3) = diff(R_Thigh_CoM(:,3))./dt; 
V_R_Thigh_CoM_Interpft(:,3) = (interpft(V_R_Thigh_CoM_Raw(:,3),nl));
V_R_Thigh_CoM(:,3) = V_R_Thigh_CoM_Interpft(:,3) - V_CoM(:,3);

    % Arm Right
V_R_Arm_CoM_Raw(:,3) = diff(R_Arm_CoM(:,3))./dt; 
V_R_Arm_CoM_Interpft(:,3) = (interpft(V_R_Arm_CoM_Raw(:,3),nl));
V_R_Arm_CoM(:,3) = V_R_Arm_CoM_Interpft(:,3) - V_CoM(:,3);

    % Forearm Right
V_R_Forearm_CoM_Raw(:,3) = diff(R_Forearm_CoM(:,3))./dt; 
V_R_Forearm_CoM_Interpft(:,3) = (interpft(V_R_Forearm_CoM_Raw(:,3),nl));
V_R_Forearm_CoM(:,3) = V_R_Forearm_CoM_Interpft(:,3) - V_CoM(:,3);


    % Foot Left
V_L_Foot_CoM_Raw(:,3) = diff(L_Foot_CoM(:,3))./dt;                      
V_L_Foot_CoM_Interpft(:,3) = (interpft(V_L_Foot_CoM_Raw(:,3),nl));       
V_L_Foot_CoM(:,3) = V_L_Foot_CoM_Interpft(:,3) - V_CoM(:,3);        

    % Shank Left
V_L_Shank_CoM_Raw(:,3) = diff(L_Shank_CoM(:,3))./dt; 
V_L_Shank_CoM_Interpft(:,3) = (interpft(V_L_Shank_CoM_Raw(:,3),nl));
V_L_Shank_CoM(:,3) = V_L_Shank_CoM_Interpft(:,3) - V_CoM(:,3);

    % Thigh Left
V_L_Thigh_CoM_Raw(:,3) = diff(L_Thigh_CoM(:,3))./dt; 
V_L_Thigh_CoM_Interpft(:,3) = (interpft(V_L_Thigh_CoM_Raw(:,3),nl));
V_L_Thigh_CoM(:,3) = V_L_Thigh_CoM_Interpft(:,3) - V_CoM(:,3);

    % Arm Left
V_L_Arm_CoM_Raw(:,3) = diff(L_Arm_CoM(:,3))./dt; 
V_L_Arm_CoM_Interpft(:,3) = (interpft(V_L_Arm_CoM_Raw(:,3),nl));
V_L_Arm_CoM(:,3) = V_L_Arm_CoM_Interpft(:,3) - V_CoM(:,3);

    % Forearm Left
V_L_Forearm_CoM_Raw(:,3) = diff(L_Forearm_CoM(:,3))./dt; 
V_L_Forearm_CoM_Interpft(:,3) = (interpft(V_L_Forearm_CoM_Raw(:,3),nl));
V_L_Forearm_CoM(:,3) = V_L_Forearm_CoM_Interpft(:,3) - V_CoM(:,3);
    
   % Medium Trunk
V_Trunk_Medium_CoM_Raw(:,3) = diff(Medium_Trunk_CoM(:,3))./dt;  
V_Trunk_Medium_CoM_Interpft(:,3) = (interpft(V_Trunk_Medium_CoM_Raw(:,3),nl));
V_Trunk_Medium_CoM(:,3) = V_Trunk_Medium_CoM_Interpft(:,3) - V_CoM(:,3);


    %Gyration Radium (GyR)
    
GyR_R_Foot = 0.690*(R_Toe - R_Ankle);            % Radius of Gyration of R_Foot
GyR_R_Shank = 0.528*(R_Ankle - R_Knee);          % Radius of Gyration of R_Shank
GyR_R_Thigh = 0.540*(R_Knee - R_Throc);          % Radius of Gyration of R_Thigh
GyR_R_Arm = 0.542*(R_Elbow - R_Shoulder);        % Radius of Gyration of R_Arm
GyR_R_Forearm = 0.526*(R_Wrist - R_Elbow);       % Radius of Gyration of R_Forearm and Hand

GyR_L_Foot = 0.690*(L_Toe - L_Ankle);            % Radius of Gyration of L_Foot
GyR_L_Shank = 0.528*(L_Ankle - L_Knee);          % Radius of Gyration of L_Shank
GyR_L_Thigh = 0.540*(L_Knee - L_Throc);          % Radius of Gyration of L_Thigh
GyR_L_Arm = 0.542*(L_Elbow - L_Shoulder);        % Radius of Gyration of L_Arm
GyR_L_Forearm = 0.526*(L_Wrist - L_Elbow);       % Radius of Gyration of L_Forearm and Hand
    
GyR_Medium_Trunk = 0.830*(Medium_Throc - Medium_Head);  % Radius of Gyration of Medium_Trunk

       
        
        %KINETIC ENERGY PER SEGMENT (LINEAR AND ANGULAR)
R_Foot_K = (0.5 .* Foot_Mass .* (V_R_Foot_CoM.^2)) + (0.5 .* Foot_Mass .* (GyR_R_Foot.^2) .* (R_Foot_Angular_Speed.^2));
R_Shank_K = (0.5 * Shank_Mass * (V_R_Shank_CoM.^2)) + (0.5 .* Shank_Mass * (GyR_R_Shank.^2) .* (R_Shank_Angular_Speed.^2));
R_Thigh_K = (0.5 * Thigh_Mass * (V_R_Thigh_CoM.^2)) + (0.5 * Thigh_Mass * (GyR_R_Thigh.^2) .* (R_Thigh_Angular_Speed.^2));
R_Arm_K = (0.5 * Arm_Mass * (V_R_Arm_CoM.^2)) + (0.5 * Arm_Mass * (GyR_R_Arm.^2) .* (R_Arm_Angular_Speed.^2));
R_Forearm_K = (0.5 * Forearm_Hand_Mass * (V_R_Forearm_CoM.^2)) + (0.5 * Forearm_Hand_Mass * (GyR_R_Forearm.^2) .* (R_Forearm_Angular_Speed.^2));

L_Foot_K = (0.5 * Foot_Mass * (V_L_Foot_CoM.^2)) + (0.5 * Foot_Mass * (GyR_L_Foot.^2) .* (L_Foot_Angular_Speed.^2));
L_Shank_K = (0.5 * Shank_Mass * (V_L_Shank_CoM.^2)) + (0.5 * Shank_Mass * (GyR_L_Shank.^2) .* (L_Shank_Angular_Speed.^2));
L_Thigh_K = (0.5 * Thigh_Mass * (V_L_Thigh_CoM.^2)) + (0.5 * Thigh_Mass * (GyR_L_Thigh.^2) .* (L_Thigh_Angular_Speed.^2));
L_Arm_K = (0.5 * Arm_Mass * (V_L_Arm_CoM.^2)) + (0.5 * Arm_Mass * (GyR_L_Arm.^2) .* (L_Arm_Angular_Speed.^2));
L_Forearm_K = (0.5 * Forearm_Hand_Mass * (V_L_Forearm_CoM.^2)) + (0.5 * Forearm_Hand_Mass * (GyR_L_Forearm.^2) .* (L_Forearm_Angular_Speed.^2));

Trunk_Medium_K = (0.5 * Trunk_Mass * (V_Trunk_Medium_CoM.^2)) + (0.5 * Trunk_Mass * (GyR_Medium_Trunk.^2) .* (Medium_Trunk_Angular_Speed.^2));



        %INTERNAL KINETIC ENERGY WITHOUT SEGMENTS TRANSFER
EK_Int = R_Foot_K + R_Shank_K + R_Thigh_K + R_Arm_K + R_Forearm_K + L_Foot_K + L_Shank_K + L_Thigh_K + L_Arm_K + L_Forearm_K + Trunk_Medium_K;
EK_Int = interpft(EK_Int, nl);
EK_Int = EK_Int(:,1) + EK_Int(:,2) + EK_Int(:,3);



        %INTERNAL KINETIC ENERGY WITH SEGMENTS TRANSFER
        
% LOWER LIMB R        
EK_Int_Limb_Lower_R = (R_Foot_K + R_Shank_K + R_Thigh_K);       % This is creating the Lower Limb R Segment
W_Int_Limb_Lower_R = diff(EK_Int_Limb_Lower_R);
W_Int_Limb_Lower_R = interpft(W_Int_Limb_Lower_R,nl);

    %STRIDE 01
W_Int_Limb_Lower_R_St1 = W_Int_Limb_Lower_R(St1_TD1:St1_TD2,:);

for i=1:(St1_TD2 - St1_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Lower_R_St1(i,j) > 0
     W_Int_Limb_Lower_R_St1_Positive(i,j) = W_Int_Limb_Lower_R_St1(i,j);    
     else
     W_Int_Limb_Lower_R_St1_Negative(i,j) = W_Int_Limb_Lower_R_St1(i,j);    
     end
   end
end

W_Int_Limb_Lower_R_St1_Positive = sum(W_Int_Limb_Lower_R_St1_Positive);

     %STRIDE 02
W_Int_Limb_Lower_R_St2 = W_Int_Limb_Lower_R(St2_TD1:St2_TD2,:);

for i=1:(St2_TD2 - St2_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Lower_R_St2(i,j) > 0
     W_Int_Limb_Lower_R_St2_Positive(i,j) = W_Int_Limb_Lower_R_St2(i,j);    
     else
     W_Int_Limb_Lower_R_St2_Negative(i,j) = W_Int_Limb_Lower_R_St2(i,j);    
     end
   end
end

W_Int_Limb_Lower_R_St2_Positive = sum(W_Int_Limb_Lower_R_St2_Positive);
 
     %STRIDE 03
W_Int_Limb_Lower_R_St3 = W_Int_Limb_Lower_R(St3_TD1:St3_TD2,:);

for i=1:(St3_TD2 - St3_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Lower_R_St3(i,j) > 0
     W_Int_Limb_Lower_R_St3_Positive(i,j) = W_Int_Limb_Lower_R_St3(i,j);    
     else
     W_Int_Limb_Lower_R_St3_Negative(i,j) = W_Int_Limb_Lower_R_St3(i,j);    
     end
   end
end

W_Int_Limb_Lower_R_St3_Positive = sum(W_Int_Limb_Lower_R_St3_Positive);


% UPPER LIMB R
EK_Int_Limb_Upper_R = (R_Arm_K + R_Forearm_K);                % This is creating the Upper Limb R Segment
W_Int_Limb_Upper_R = diff(EK_Int_Limb_Upper_R);
W_Int_Limb_Upper_R = interpft(W_Int_Limb_Upper_R,nl);

     %STRIDE 01
W_Int_Limb_Upper_R_St1 = W_Int_Limb_Upper_R(St1_TD1:St1_TD2,:);

for i=1:(St1_TD2 - St1_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Upper_R_St1(i,j) > 0
     W_Int_Limb_Upper_R_St1_Positive(i,j) = W_Int_Limb_Upper_R_St1(i,j);    
     else
     W_Int_Limb_Upper_R_St1_Negative(i,j) = W_Int_Limb_Upper_R_St1(i,j);    
     end
   end
end
 W_Int_Limb_Upper_R_St1_Positive = sum(W_Int_Limb_Upper_R_St1_Positive);
 
    %STRIDE 02
W_Int_Limb_Upper_R_St2 = W_Int_Limb_Upper_R(St2_TD1:St2_TD2,:);

for i=1:(St2_TD2 - St2_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Upper_R_St2(i,j) > 0
     W_Int_Limb_Upper_R_St2_Positive(i,j) = W_Int_Limb_Upper_R_St2(i,j);    
     else
     W_Int_Limb_Upper_R_St2_Negative(i,j) = W_Int_Limb_Upper_R_St2(i,j);    
     end
   end
end
 W_Int_Limb_Upper_R_St2_Positive = sum(W_Int_Limb_Upper_R_St2_Positive);

     %STRIDE 03
W_Int_Limb_Upper_R_St3 = W_Int_Limb_Upper_R(St3_TD1:St3_TD2,:);

for i=1:(St3_TD2 - St3_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Upper_R_St3(i,j) > 0
     W_Int_Limb_Upper_R_St3_Positive(i,j) = W_Int_Limb_Upper_R_St3(i,j);    
     else
     W_Int_Limb_Upper_R_St3_Negative(i,j) = W_Int_Limb_Upper_R_St3(i,j);    
     end
   end
end
 W_Int_Limb_Upper_R_St3_Positive = sum(W_Int_Limb_Upper_R_St3_Positive);
 

% LOWER LIMB L
EK_Int_Limb_Lower_L = (L_Foot_K + L_Shank_K + L_Thigh_K);       % This is creating the Lower Limb L Segment
W_Int_Limb_Lower_L = diff(EK_Int_Limb_Lower_L);
W_Int_Limb_Lower_L = interpft(W_Int_Limb_Lower_L,nl);


    %STRIDE 01
W_Int_Limb_Lower_L_St1 = W_Int_Limb_Lower_L(St1_TD1:St1_TD2,:);
    
for i=1:(St1_TD2 - St1_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Lower_L_St1(i,j) > 0
     W_Int_Limb_Lower_L_St1_Positive(i,j) = W_Int_Limb_Lower_L_St1(i,j);    
     else
     W_Int_Limb_Lower_L_St1_Negative(i,j) = W_Int_Limb_Lower_L_St1(i,j);    
     end
   end
end
 W_Int_Limb_Lower_L_St1_Positive = sum(W_Int_Limb_Lower_L_St1_Positive);

    %STRIDE 02
W_Int_Limb_Lower_L_St2 = W_Int_Limb_Lower_L(St2_TD1:St2_TD2,:);
    
for i=1:(St2_TD2 - St2_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Lower_L_St2(i,j) > 0
     W_Int_Limb_Lower_L_St2_Positive(i,j) = W_Int_Limb_Lower_L_St2(i,j);    
     else
     W_Int_Limb_Lower_L_St2_Negative(i,j) = W_Int_Limb_Lower_L_St2(i,j);    
     end
   end
end
 W_Int_Limb_Lower_L_St2_Positive = sum(W_Int_Limb_Lower_L_St2_Positive);

   %STRIDE 03
W_Int_Limb_Lower_L_St3 = W_Int_Limb_Lower_L(St3_TD1:St3_TD2,:);
    
for i=1:(St3_TD2 - St3_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Lower_L_St3(i,j) > 0
     W_Int_Limb_Lower_L_St3_Positive(i,j) = W_Int_Limb_Lower_L_St3(i,j);    
     else
     W_Int_Limb_Lower_L_St3_Negative(i,j) = W_Int_Limb_Lower_L_St3(i,j);    
     end
   end
end
 W_Int_Limb_Lower_L_St3_Positive = sum(W_Int_Limb_Lower_L_St3_Positive);

 
% UPPER LIMB L
EK_Int_Limb_Upper_L = (L_Arm_K + L_Forearm_K);                % This is creating the Upper Limb L Segment
W_Int_Limb_Upper_L = diff(EK_Int_Limb_Upper_L);
W_Int_Limb_Upper_L = interpft(W_Int_Limb_Upper_L,nl);


  %STRIDE 01
W_Int_Limb_Upper_L_St1 = W_Int_Limb_Upper_L(St1_TD1:St1_TD2,:);

for i=1:(St1_TD2 - St1_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Upper_L_St1(i,j) > 0
     W_Int_Limb_Upper_L_St1_Positive(i,j) = W_Int_Limb_Upper_L_St1(i,j);    
     else
     W_Int_Limb_Upper_L_St1_Negative(i,j) = W_Int_Limb_Upper_L_St1(i,j);    
     end
   end
end
 W_Int_Limb_Upper_L_St1_Positive = sum(W_Int_Limb_Upper_L_St1_Positive);

  %STRIDE 02
W_Int_Limb_Upper_L_St2 = W_Int_Limb_Upper_L(St2_TD1:St2_TD2,:);

for i=1:(St2_TD2 - St2_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Upper_L_St2(i,j) > 0
     W_Int_Limb_Upper_L_St2_Positive(i,j) = W_Int_Limb_Upper_L_St2(i,j);    
     else
     W_Int_Limb_Upper_L_St2_Negative(i,j) = W_Int_Limb_Upper_L_St2(i,j);    
     end
   end
end
 W_Int_Limb_Upper_L_St2_Positive = sum(W_Int_Limb_Upper_L_St2_Positive);

   %STRIDE 03
W_Int_Limb_Upper_L_St3 = W_Int_Limb_Upper_L(St3_TD1:St3_TD2,:);

for i=1:(St3_TD2 - St3_TD1 + 1)
   for j=1:3
     if W_Int_Limb_Upper_L_St3(i,j) > 0
     W_Int_Limb_Upper_L_St3_Positive(i,j) = W_Int_Limb_Upper_L_St3(i,j);    
     else
     W_Int_Limb_Upper_L_St3_Negative(i,j) = W_Int_Limb_Upper_L_St3(i,j);    
     end
   end
end
 W_Int_Limb_Upper_L_St3_Positive = sum(W_Int_Limb_Upper_L_St3_Positive);

 
% TRUNK
EK_Int_Trunk = Trunk_Medium_K;                               % This is creating the Trunk Segment
W_Int_Trunk = diff(EK_Int_Trunk);
W_Int_Trunk = interpft(W_Int_Trunk,nl);


    %STRIDE 01
W_Int_Trunk_St1 = W_Int_Trunk(St1_TD1:St1_TD2,:);

for i=1:(St1_TD2 - St1_TD1 + 1)
   for j=1:3
     if W_Int_Trunk_St1(i,j) > 0
     W_Int_Trunk_St1_Positive(i,j) = W_Int_Trunk_St1(i,j);    
     else
     W_Int_Trunk_St1_Negative(i,j) = W_Int_Trunk_St1(i,j);    
     end
   end
end
 W_Int_Trunk_St1_Positive = sum(W_Int_Trunk_St1_Positive);

    %STRIDE 02
W_Int_Trunk_St2 = W_Int_Trunk(St2_TD1:St2_TD2,:);

for i=1:(St2_TD2 - St2_TD1 + 1)
   for j=1:3
     if W_Int_Trunk_St2(i,j) > 0
     W_Int_Trunk_St2_Positive(i,j) = W_Int_Trunk_St2(i,j);    
     else
     W_Int_Trunk_St2_Negative(i,j) = W_Int_Trunk_St2(i,j);    
     end
   end
end
 W_Int_Trunk_St2_Positive = sum(W_Int_Trunk_St2_Positive);

    %STRIDE 03
W_Int_Trunk_St3 = W_Int_Trunk(St3_TD1:St3_TD2,:);

for i=1:(St3_TD2 - St3_TD1 + 1)
   for j=1:3
     if W_Int_Trunk_St3(i,j) > 0
     W_Int_Trunk_St3_Positive(i,j) = W_Int_Trunk_St3(i,j);    
     else
     W_Int_Trunk_St3_Negative(i,j) = W_Int_Trunk_St3(i,j);    
     end
   end
end
 W_Int_Trunk_St3_Positive = sum(W_Int_Trunk_St3_Positive);
 
 
%% SECTION 07
% TOTAL MECHANICAL ENERGY AND WORK %
    
EM = Pg + EK_Total + EK_Int;
EM_Ext = Pg + EK_Total;

EM_St1 = EM(St1_TD1:St1_TD2,:);
EM_St2 = EM(St2_TD1:St2_TD2,:);
EM_St3 = EM(St3_TD1:St3_TD2,:);

EM_Ext_St1 = EM_Ext(St1_TD1:St1_TD2,:);
EM_Ext_St2 = EM_Ext(St2_TD1:St2_TD2,:);
EM_Ext_St3 = EM_Ext(St3_TD1:St3_TD2,:);

 
W_Int_Transfer_St1 = W_Int_Limb_Lower_R_St1_Positive + W_Int_Limb_Upper_R_St1_Positive + W_Int_Limb_Lower_L_St1_Positive + W_Int_Limb_Upper_L_St1_Positive + W_Int_Trunk_St1_Positive;
W_Int_Transfer_St1 = sum(W_Int_Transfer_St1);
W_Int_Transfer_St2 = W_Int_Limb_Lower_R_St2_Positive + W_Int_Limb_Upper_R_St2_Positive + W_Int_Limb_Lower_L_St2_Positive + W_Int_Limb_Upper_L_St2_Positive + W_Int_Trunk_St2_Positive;
W_Int_Transfer_St2 = sum(W_Int_Transfer_St2);
W_Int_Transfer_St3 = W_Int_Limb_Lower_R_St3_Positive + W_Int_Limb_Upper_R_St3_Positive + W_Int_Limb_Lower_L_St3_Positive + W_Int_Limb_Upper_L_St3_Positive + W_Int_Trunk_St3_Positive;
W_Int_Transfer_St3 = sum(W_Int_Transfer_St3);


    %MECHANICAL WORK (J) AND MECHANICAL POWER(W)

    
    %TOTAL MECHANICAL WORK AND POWER
W_Total = diff(EM);
W_Total = sum(W_Total(W_Total>0));
Power_Mechanical_Total = W_Total/(length(Kinematic)*dt);

W_Ext_Total = diff(EM_Ext);
W_Ext_Total_Positive = sum(W_Ext_Total(W_Ext_Total>0));
Power_Ext_Mechanical_Total = W_Ext_Total_Positive/(length(Kinematic)*dt);


    %Stride 01
W_St1 = diff(EM_St1);
W_St1 = sum(W_St1(W_St1>0));    
Power_Mechanical_St1 = W_St1/(length(EM_St1)*dt);

W_Ext_St1 = diff(EM_Ext_St1);
W_Ext_St1 = sum(W_Ext_St1(W_Ext_St1>0));    
Power_Ext_Mechanical_St1 = W_Ext_St1/(length(EM_Ext_St1)*dt);

W_Transfer_St1 = W_Ext_St1 + W_Int_Transfer_St1;
Power_Mechanical_Transfer_St1 = W_Transfer_St1/(length(EM_St1)*dt);

    %Stride 02
W_St2 = diff(EM_St2);
W_St2 = sum(W_St2(W_St2>0));  
Power_Mechanical_St2 = W_St2/(length(EM_St2)*dt);


W_Ext_St2 = diff(EM_Ext_St2);
W_Ext_St2 = sum(W_Ext_St2(W_Ext_St2>0));  
Power_Ext_Mechanical_St2 = W_Ext_St2/(length(EM_Ext_St2)*dt);

W_Transfer_St2 = W_Ext_St2 + W_Int_Transfer_St2;
Power_Mechanical_Transfer_St2 = W_Transfer_St2/(length(EM_St2)*dt);

    %Stride 03
W_St3 = diff(EM_St3);
W_St3 = sum(W_St3(W_St3>0));  
Power_Mechanical_St3 = W_St3/(length(EM_St3)*dt);

W_Ext_St3 = diff(EM_Ext_St3);
W_Ext_St3 = sum(W_Ext_St3(W_Ext_St3>0));  
Power_Ext_Mechanical_St3 = W_Ext_St3/(length(EM_Ext_St3)*dt);

W_Transfer_St3 = W_Ext_St3 + W_Int_Transfer_St3;
Power_Mechanical_Transfer_St3 = W_Transfer_St3/(length(EM_St3)*dt);

%% SECTION 08
% EXPORT DATA 

Output_File_Path_Full = [Output_File_Path '\' File_Name '_EM' '.xls'];

Export_EM = [EM_Ext];
Export_Header_EM = {'EM_Ext'};
Export_Excel_EM = [Export_Header_EM; num2cell(Export_EM)];
xlswrite (Output_File_Path_Full, Export_Excel_EM);             

Export_Variables = [W_St1 W_St2 W_St3 Power_Mechanical_St1 Power_Mechanical_St2 Power_Mechanical_St3 W_Ext_St1 W_Ext_St2 W_Ext_St3 Power_Ext_Mechanical_St1 Power_Ext_Mechanical_St2 Power_Ext_Mechanical_St3 W_Transfer_St1 W_Transfer_St2 W_Transfer_St3 Power_Mechanical_Transfer_St1 Power_Mechanical_Transfer_St2 Power_Mechanical_Transfer_St3 St1_Duration St2_Duration St3_Duration St1_Frequency St2_Frequency St3_Frequency St1_Lenght St2_Lenght St3_Lenght St1_Speed St2_Speed St3_Speed St1_Contact_Phase_Duration St2_Contact_Phase_Duration St3_Contact_Phase_Duration St1_Duty_Factor St2_Duty_Factor St3_Duty_Factor];
Export_Header = {'W_St1', 'W_St2', 'W_St3', 'Power_Mechanical_St1', 'Power_Mechanical_St2', 'Power_Mechanical_St3', 'W_Ext_St1', 'W_Ext_St2', 'W_Ext_St3', 'Power_Ext_Mechanical_St1', 'Power_Ext_Mechanical_St2', 'Power_Ext_Mechanical_St3', 'W_Transfer_St1', 'W_Transfer_St2', 'W_Transfer_St3', 'Power_Mechanical_Transfer_St1', 'Power_Mechanical_Transfer_St2', 'Power_Mechanical_Transfer_St3', 'St1_Duration', 'St2_Duration', 'St3_Duration', 'St1_Frequency', 'St2_Frequency', 'St3_Frequency', 'St1_Lenght', 'St2_Lenght', 'St3_Lenght', 'St1_Speed', 'St2_Speed', 'St3_Speed', 'St1_Contact_Phase_Duration', 'St2_Contact_Phase_Duration', 'St3_Contact_Phase_Duration', 'St1_Duty_Factor', 'St2_Duty_Factor', 'St3_Duty_Factor'};
Export_Excel = [Export_Header; num2cell(Export_Variables)];
xlswrite (Output_File_Path_Full, Export_Excel);             

%% SECTION 09
% GRAPHICS 

figure ('Name','SHANK POSITIONS')
title ('SHANK  Positions');
xlabel('Frames (n)');
ylabel ('Position (m)');
hold on
plot(R_Shank(:,1),'r');
hold on
plot(R_Shank(:,2),'b');
hold on
plot(R_Shank(:,3),'g');
hold on
legend('Shank X','Shank Y','Shank Z')


figure ('Name','CoM Position Y')
title ('CoM Position Y');
xlabel('Frames (n)');
ylabel ('Position (m)');
hold on
plot (CoM(:,2),'g:','LineWidth',3.0);
hold on
legend ('CoM Y Position')


figure ('Name','Segment CoM')
title ('Segment CoM');
xlabel('Frames (n)');
ylabel ('Position (m)');
hold on
plot (R_Thigh_CoM(:,2),'k:','LineWidth',3.0);
hold on
plot (L_Thigh_CoM(:,2),'r:','LineWidth',3.0);
hold on
hold on
plot (R_Shank_CoM(:,2),'b:','LineWidth',3.0);
hold on
plot (L_Shank_CoM(:,2),'g:','LineWidth',3.0);
legend ('R Thigh','L Thigh','R Shank','L Shank')

figure ('Name','Potential Gravitational Energy')
title ('Potential Gravitational Energy');
xlabel('Frames (n)');
ylabel ('Potential Gravitational Energy (J)');
hold on
plot (Pg,'r');


figure ('Name','Kinetic Energy Total')
title ('Kinetic Energy');
xlabel('Frames (n)');
ylabel ('Kinetic Energy (J)');
hold on
plot (EK_Total,'r:');
hold on
legend ('EK')


figure ('Name','Kinetic Energy X')
title ('Kinetic Energy');
xlabel('Frames (n)');
ylabel ('Kinetic Energy (J)');
hold on
plot (EK_X,'r');
hold on
legend ('EK_X')


figure ('Name','Kinetic Energy Y')
title ('Kinetic Energy');
xlabel('Frames (n)');
ylabel ('Kinetic Energy (J)');
plot (EK_Y,'b');
hold on
legend ('EK_Y')


figure ('Name','Kinetic Energy Z')
title ('Kinetic Energy');
xlabel('Frames (n)');
ylabel ('Kinetic Energy (J)');
hold on
plot (EK_Z,'g');
hold on
legend ('EK_Z')


figure('Name','Internal Kinetic Energy')
title ('Internal Kinetic Energy');
xlabel('Frames (n)');
ylabel ('Energy (J)');
hold on
plot(EK_Int,'r');


figure ('Name','External Mechanical Energy')
title ('External Mechanical Energy');
xlabel('Time (s)');
ylabel ('Mechanical Energy (J)');
hold on
plot (Time_Row(St1_TD1:St1_TD2), EM_Ext(St1_TD1:St1_TD2), 'r', 'LineWidth', 4);
hold on
plot (Time_Row(St2_TD1:St2_TD2), EM_Ext(St2_TD1:St2_TD2), 'g', 'LineWidth', 4);
hold on
plot (Time_Row(St3_TD1:St3_TD2), EM_Ext(St3_TD1:St3_TD2), 'b', 'LineWidth', 4);
legend ('EM')


figure ('Name','Total Mechanical Energy')
title ('Total Mechanical Energy');
xlabel('Time (s)');
ylabel ('Mechanical Energy (J)');
hold on
plot (Time_Row(St1_TD1:St1_TD2), EM(St1_TD1:St1_TD2), 'r', 'LineWidth', 4);
hold on
plot (Time_Row(St2_TD1:St2_TD2), EM(St2_TD1:St2_TD2), 'g', 'LineWidth', 4);
hold on
plot (Time_Row(St3_TD1:St3_TD2), EM(St3_TD1:St3_TD2), 'b', 'LineWidth', 4);


y = [Power_Mechanical_St1 Power_Mechanical_Transfer_St1; Power_Mechanical_St2 Power_Mechanical_Transfer_St2; Power_Mechanical_St3 Power_Mechanical_Transfer_St3];
figure ('Name','Mechanical Power')
title ('Mechanical Power Without and With Segments Transfer');
xlabel('Stride');
ylabel ('Mechanical Power (W)')
hold on
bar ([1 2 3], y)
legend('Power Without', 'Power With')