%データファイルから二次元プロットするサンプル
function plot_2dgraph_datafile
    %初期化
    close all;
    clear;
    %グラフの表示範囲指定[下限値：刻み幅：上限値]
    xmin = 0;
    xmax = 200;
    ymin = 0;
    ymax = 280; 
    graph_meshwidth_x = 10.0;
    graph_meshwidth_y = 10.0;
    
    %-----------------------------------------------------
    %グラフ描画用データ作成
    %----------------------------------------------------- 

    M = dlmread('result.csv');%データを読み込む
    %以後ここに列の数を増やすことができる

    %-----------------------------------------------------
    %グラフ描画
    %-----------------------------------------------------
    t = linspace(0,2*pi,100);
    hold on;
    for I = transpose(M)
        cx = I(1); 
        cy = I(2); % 中心
        r = I(3);           % 半径
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
    set(gca, 'xtick',xmin:graph_meshwidth_x:xmax); %メモリの刻み方
    set(gca, 'ytick',ymin:graph_meshwidth_y:ymax); %メモリの刻み方
    ylabel('x,yの値','FontSize',20);
    xlabel('tの値','FontSize',20);
    
    saveas(gcf,'output.eps' , 'epsc2'); %epsの色出力
    
end

