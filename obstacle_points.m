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

    M = dlmread('env_simple.txt');%�f�[�^��ǂݍ���
    X =M(:,1); %1��ڂ�ǂݍ���
    Y =M(:,2); %2��ڂ�ǂݍ���
    %�Ȍケ���ɗ�̐��𑝂₷���Ƃ��ł���
    
    %-----------------------------------------------------
    %�O���t�`��
    %-----------------------------------------------------
    
    scatter(X,Y);%plot(T,X1,T,X2,T,X3,T,X4)�Ȃǂ̂悤�ɏo�͂���O���t�̐��𑝂₷���Ƃ��ł���
    
    axis([xmin, xmax, ymin, ymax]);
    set(gca, 'xtick',xmin:graph_meshwidth_x:xmax); %�������̍��ݕ�
    set(gca, 'ytick',ymin:graph_meshwidth_y:ymax); %�������̍��ݕ�
    ylabel('Y�̒l','FontSize',20);
    xlabel('X�̒l','FontSize',20);
    %legend('x = 1+t','y = 2+t','Location','Best')%�}��̐ݒ�,'Best'��'SouthWest'��'NorthEast'�ȂǂɕύX�\
    
    saveas(gcf,'output.eps' , 'epsc2'); %eps�̐F�o��
    
end

