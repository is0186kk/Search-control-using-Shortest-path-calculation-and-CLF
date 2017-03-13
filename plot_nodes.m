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
    graph_meshwidth_x = 10.0;
    graph_meshwidth_y = 10.0;
    
    %-----------------------------------------------------
    %�O���t�`��p�f�[�^�쐬
    %----------------------------------------------------- 

    M = dlmread('result.csv');%�f�[�^��ǂݍ���
    %�Ȍケ���ɗ�̐��𑝂₷���Ƃ��ł���

    %-----------------------------------------------------
    %�O���t�`��
    %-----------------------------------------------------
    t = linspace(0,2*pi,100);
    hold on;
    for I = transpose(M)
        cx = I(1); 
        cy = I(2); % ���S
        r = I(3);           % ���a
        n = I(4)+1;
        switch n
        case 1
            plot(r*sin(t)+cx,r*cos(t)+cy,'k');
        case 2
            plot(r*sin(t)+cx,r*cos(t)+cy,'b');
        case 3
            plot(r*sin(t)+cx,r*cos(t)+cy,'g');
        case 3
            plot(r*sin(t)+cx,r*cos(t)+cy,'r');
        end
    end
    hold off;
    
    axis([xmin, xmax, ymin, ymax]);
    set(gca, 'xtick',xmin:graph_meshwidth_x:xmax); %�������̍��ݕ�
    set(gca, 'ytick',ymin:graph_meshwidth_y:ymax); %�������̍��ݕ�
    ylabel('x,y�̒l','FontSize',20);
    xlabel('t�̒l','FontSize',20);
    
    saveas(gcf,'output.eps' , 'epsc2'); %eps�̐F�o��
    
end

