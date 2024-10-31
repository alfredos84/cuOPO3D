clear all; close all; clc;


N=0.5;
NX = 64; NY = NX;
LX = 1e-3; LY = LX;
dx = LX/(NX-1); dy = dx;

x = linspace(-LX/2, LX/2, NX);
y = linspace(-LY/2, LY/2, NY);

for s = 0

    Ar = load (['test/sPPLT_N_',num2str(N),'/inic_slice_XY_p_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/inic_slice_XY_p_',num2str(s),'_i.dat']);


    A1 = Ar + 1i*Ai;

    subplot(2,3,1)
    imagesc(x,y,abs(A1).^2)
    % zlim   ([0,6])
    % shading interp;
    title('Input pump')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];

    %% M2
    Ar = load (['test/sPPLT_N_',num2str(N),'/M2_s_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/M2_s_',num2str(s),'_i.dat']);

    M2 = Ar + 1i*Ai;

    subplot(2,3,2)
    imagesc(x,y,abs(M2).^2)
    % zlim   ([0,6])
    % shading interp;
    title('M2')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];


    %% M3
    Ar = load (['test/sPPLT_N_',num2str(N),'/M3_s_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/M3_s_',num2str(s),'_i.dat']);

    M3 = Ar + 1i*Ai;

    subplot(2,3,3)
    imagesc(x,y,abs(M3).^2)
    % zlim   ([0,6])
    % shading interp;
    title('M3')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];


    %% M1_s_
    Ar = load (['test/sPPLT_N_',num2str(N),'/M1_s_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/M1_s_',num2str(s),'_i.dat']);

    M1 = Ar + 1i*Ai;

    subplot(2,3,4)
    imagesc(x,y,abs(M1).^2)
    % zlim   ([0,6])
    % shading interp;
    title('M3')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];


    %% InputAgain
    Ar = load (['test/sPPLT_N_',num2str(N),'/InputAgain_s_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/InputAgain_s_',num2str(s),'_i.dat']);

    In = Ar + 1i*Ai;

    subplot(2,3,5)
    imagesc(x,y,abs(In).^2)
    % zlim   ([0,6])
    % shading interp;
    title('Input Again')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];

    subplot(2,3,6)
    imagesc( x, y, abs(In).^2 - abs(A1).^2 )
    % zlim   ([0,6])
    % shading interp;
    title('Comparison')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];

    % waitforbuttonpress
end


figure();

for s = 0

    Ar = load (['test/sPPLT_N_',num2str(N),'/inic_slice_XY_p_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/inic_slice_XY_p_',num2str(s),'_i.dat']);


    A1 = Ar + 1i*Ai;

    subplot(2,3,1)
    imagesc(x,y,angle(A1)/pi)
    % zlim   ([0,6])
    % shading interp;
    title('Input pump')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];

    %% M2
    Ar = load (['test/sPPLT_N_',num2str(N),'/M2_s_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/M2_s_',num2str(s),'_i.dat']);

    M2 = Ar + 1i*Ai;

    subplot(2,3,2)
    imagesc(x,y,angle(M2)/pi)
    % zlim   ([0,6])
    % shading interp;
    title('M2')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];


    %% M3
    Ar = load (['test/sPPLT_N_',num2str(N),'/M3_s_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/M3_s_',num2str(s),'_i.dat']);

    M3 = Ar + 1i*Ai;

    subplot(2,3,3)
    imagesc(x,y,angle(M3)/pi)
    % zlim   ([0,6])
    % shading interp;
    title('M3')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];


    %% M1_s_
    Ar = load (['test/sPPLT_N_',num2str(N),'/M1_s_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/M1_s_',num2str(s),'_i.dat']);

    M1 = Ar + 1i*Ai;

    subplot(2,3,4)
    imagesc(x,y,angle(M1)/pi)
    % zlim   ([0,6])
    % shading interp;
    title('M1')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];


    %% InputAgain
    Ar = load (['test/sPPLT_N_',num2str(N),'/InputAgain_s_',num2str(s),'_r.dat']);
    Ai = load (['test/sPPLT_N_',num2str(N),'/InputAgain_s_',num2str(s),'_i.dat']);

    In = Ar + 1i*Ai;

    subplot(2,3,5)
    imagesc(x,y,angle(In)/pi)
    % zlim   ([0,6])
    % shading interp;
    title('Input Again')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];

    subplot(2,3,6)
    imagesc( x, y, angle(In)/pi - angle(A1)/pi )
    % zlim   ([0,6])
    % shading interp;
    title('Comparison')
    colorbar();
    ax = gca;
    ax.PlotBoxAspectRatio = [1,1,1];

    % waitforbuttonpress
end