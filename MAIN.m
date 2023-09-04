clc,clear;

%% MAIN
% This code provide reliability analysis of HEC-RAS 3D model.
% 
% HEC-RAS model is not provided in this bundle but the connection code
% between Kriging reliability method and HEC-RAS model using MATLAB is
% provided
% 
% For Kriging model, Kriging algorithm used in this study is solely taken 
% from DACE toolbox Thanks to: Nielsen, H. B. (Author), Lophaven, S. N. 
% (Author), & Søndergaard, J. (Author). (2002). DACE - A Matlab Kriging 
% Toolbox. Computer programme, Informatics and Mathematical Modelling, 
% Technical University of Denmark, DTU. 
% http://www2.imm.dtu.dk/pubdb/p.php?1460
% 
% For Kriging Reliability method, AK-PSO-HHs is utilized 
% Adaptive Kriging Adopting PSO with Hollow-Hypersphere space in structural 
% reliability assessmen 
% Author John Thedy, Kuo-Wei Liao
% https://doi.org/10.1016/j.probengmech.2023.103513

%% AK-PSO-HHs initial parameter
% Mean, std, and random variable distribution of problem could be defined
% in problem.m

% In this example case, there are four ramdon variables which is rainfall
% intensity, manning coefficient, skewness and kurtosis of pearson pdf.
% skewness and kurtosis of pearson pdf is used to defined rainfall
% distribution along raining duration.

% dmin is minimal distance between DoE

% uxdoe and lxdoe is upper and lower bound of DoE in standard normal space

% initsample defined as initial sample size to create Kriging

% beta1 and beta2 area initial hyper hollow sphere outer and inner diameter
% respectively

probtype=1;
dmin=0.01;
uxdoe=8;lxdoe=-8;
initsample=20;
beta1=8;beta2=3;

%% HEC-RAS model input file

% In this section, user need to create and modified u file and g file from
% HEC-RAS main file. Example of modified u(Ori_ok.u33) and g(Ori_ok.g03) 
% file are provided.

% UOriginp is original input file containing rainfall information
% Unewinp is modified UOriginp to create new file with new rainfall
UOriginp='D:\OneDrive\Documents\Work\HECRAS2\230506 baolixi small no_gsw\Hec-ras Baolixi 230505 _no gsw\Ori_ok.u33';
Unewinp='D:\OneDrive\Documents\Work\HECRAS2\230506 baolixi small no_gsw\Hec-ras Baolixi 230505 _no gsw\ok.u33';

% GOriginp is original input file containing Manning information
% Gnewinp is modified UOriginp to create new file with new Manning coeff
GOriginp='D:\OneDrive\Documents\Work\HECRAS2\230506 baolixi small no_gsw\Hec-ras Baolixi 230505 _no gsw\Ori_ok.g03';
Gnewinp='D:\OneDrive\Documents\Work\HECRAS2\230506 baolixi small no_gsw\Hec-ras Baolixi 230505 _no gsw\ok.g03';

%PRJfile is HEC-RAS main file
%HDFfile is HEC-RAS output file
PRJfile='D:\OneDrive\Documents\Work\HECRAS2\230506 baolixi small no_gsw\Hec-ras Baolixi 230505 _no gsw\ok.prj';
HDFfile='D:\OneDrive\Documents\Work\HECRAS2\230506 baolixi small no_gsw\Hec-ras Baolixi 230505 _no gsw\ok.p81.hdf';

% Final Result Contain failure probability, number of sample, and MPP
Result=AKPSOHHS(probtype,dmin,uxdoe,lxdoe,initsample,beta1,beta2,UOriginp,Unewinp,GOriginp,Gnewinp,PRJfile,HDFfile);
