close all 
clear all
clc

%pin A 
syms R_Ax R_Ay R_AH R_AB
theta_1=atand(2/2.5);
Eq_01=R_Ax==R_AB+R_AH*cosd(theta_1);
Eq_02=R_Ay==R_AB*sind(theta_1);

%pin B
syms R_BH R_BI R_BC
theta_2=atand(2/2.5);
Eq_03=R_AB+R_BH*cosd(theta_2)==R_BC;
Eq_04=R_BI+R_BH*sind(theta_2)==0;

%pin C
syms R_CD R_CI R_CJ
theta_3=atand(3.5/5.5);
Eq_05=R_BC+R_CI*cosd(theta_3)==R_CD;
Eq_06=R_CJ+R_CI*sind(theta_3)==0;

%pin D
syms R_DE R_DJ R_DK
theta_4=atand(5/6);
Eq_07=R_CD+R_DJcosd(theta_4)==R_DL*cosd(theta_4)==R_DE;
Eq_08=(R_DJ+R_DL)*sind(theta_4)+R_DK+100==0;

%pin E
syms R_LE R_EM R_EF
Eq_09=R_EF+R_EM*cosd(theta_3)==R_DE;
Eq_10=R_LE+R_EM*sind(theta_3)==0;

%pin F
syms R_MF R_NF R_GF
Eq_09=R_GF+R_NF*cosd(theta_3)==R_EF;
Eq_10=R_MF+R_NF*sind(theta_3)==0;

%pin G
syms R_NG R_Gy
Eq_13=R_GF+R_NG*cosd(theta_1)==0;
Eq_14=R_Gy==R_NG*sind(theta_1);

%pin h
syms R_IH
theta_5=atand(1.5/2.5);
Eq_15= R_AH*cosd(theta_1)==R_IH*;

Eq_16=

%pin I
syms R_IJ 
theta_6=atand(1.5/5.5);
Eq_17=
Eq_18=


%pin J
syms R_JK
theta_7=atand(.05/6);
Eq_19=
Eq_20=

%pin K
syms R_KL
theta_8=atand(3/6);
Eq_21=
Eq_22=


%pin L  
syms R_LE
theta_9=atand(3.5/5.5);
Eq_23=
Eq_24=

%pin M
syms R_MN   
theta_10=atand(2/2.5);
Eq_25=
Eq_26=

%pin N
syms R_NO
theta_11=atand(2/2.5);
Eq_27=
Eq_28=


eqns= [Eq_01 Eq_02 Eq_03 Eq_04 Eq_05 Eq_06 Eq_07 Eq_08 Eq_09 Eq_10...
    Eq_11 Eq_12 Eq_13 Eq_14 Eq_15 Eq_16 Eq_17 Eq_18 Eq_19 Eq_20...
    Eq_21 Eq_22 Eq_23 Eq_24 Eq_25 Eq_26 Eq_27 Eq_28];

vars= [R_Ax R_Ay R_AH R_AB R_BH R_BI R_BC R_CD R_CI R_CJ R_DE... 
R_DJ R_DK R_DL R_LE R_EM R_EF R_MF R_NF R_GF R_NG R_Gy R_IH R_IJ R_JK R_KL R_MN R_NO];

[sol_R_Ax,sol_R_Ay,sol_R_AH,sol_R_AB,sol_R_BH,sol_R_BI,sol_R_BC,...
    sol_R_CD,sol_R_CI,sol_R_CJ,sol_R_DE,sol_R_DJ,sol_R_DK,sol_R_DL,...
    sol_R_LE,sol_R_EM,sol_R_EF,sol_R_MF,sol_R_NF,sol_R_GF,...
    sol_R_NG,sol_R_Gy,sol_R_IH,sol_R_IJ,sol_R_JK,sol_R_KL,...
    sol_R_MN,sol_R_NO]=...
    solve(eqns,vars);