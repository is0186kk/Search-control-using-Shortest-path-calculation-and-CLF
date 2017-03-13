%�f�[�^�t�@�C������񎟌��v���b�g����T���v��
function plot_2dgraph_datafile
    %������
    close all;
    clear;
    %�O���t�̕\���͈͎w��[�����l�F���ݕ��F����l]
    xmin = 0;
    xmax = 200;
    ymin = 0;
    ymax = 280; 
    graph_meshwidth_x = 2.0;
    graph_meshwidth_y = 5.0;
    
    %-----------------------------------------------------
    %�O���t�`��p�f�[�^�쐬
    %----------------------------------------------------- 

    M = dlmread('resultXY.csv');%�f�[�^��ǂݍ���
    X =M(:,1); %1��ڂ�ǂݍ���
    Y =M(:,2); %2��ڂ�ǂݍ���
    T =M(:,3); %3��ڂ�ǂݍ���
    %�Ȍケ���ɗ�̐��𑝂₷���Ƃ��ł���
    M2 = dlmread('env_simple.txt');
    OX =M2(:,1); %1��ڂ�ǂݍ���
    OY =M2(:,2); %2��ڂ�ǂݍ���
    %-----------------------------------------------------
    %�O���t�`��
    %-----------------------------------------------------
    hold on;
    plot(X,Y);%plot(T,X1,T,X2,T,X3,T,X4)�Ȃǂ̂悤�ɏo�͂���O���t�̐��𑝂₷���Ƃ��ł���
    scatter(OX,OY);
    hold off;
    axis([xmin, xmax, ymin, ymax]);
    set(gca, 'xtick',xmin:graph_meshwidth_x:xmax); %�������̍��ݕ�
    set(gca, 'ytick',ymin:graph_meshwidth_y:ymax); %�������̍��ݕ�
    ylabel('x,y�̒l','FontSize',20);
    xlabel('t�̒l','FontSize',20);
    legend('x = 1+t','y = 2+t','Location','Best')%�}��̐ݒ�,'Best'��'SouthWest'��'NorthEast'�ȂǂɕύX�\
    
    saveas(gcf,'output.eps' , 'epsc2'); %eps�̐F�o��
    
end

