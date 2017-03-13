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
    graph_meshwidth_x = 2.0;
    graph_meshwidth_y = 5.0;
    
    %-----------------------------------------------------
    %グラフ描画用データ作成
    %----------------------------------------------------- 

    M = dlmread('env_simple.txt');%データを読み込む
    X =M(:,1); %1列目を読み込む
    Y =M(:,2); %2列目を読み込む
    %以後ここに列の数を増やすことができる
    
    %-----------------------------------------------------
    %グラフ描画
    %-----------------------------------------------------
    
    scatter(X,Y);%plot(T,X1,T,X2,T,X3,T,X4)などのように出力するグラフの数を増やすこともできる
    
    axis([xmin, xmax, ymin, ymax]);
    set(gca, 'xtick',xmin:graph_meshwidth_x:xmax); %メモリの刻み方
    set(gca, 'ytick',ymin:graph_meshwidth_y:ymax); %メモリの刻み方
    ylabel('Yの値','FontSize',20);
    xlabel('Xの値','FontSize',20);
    %legend('x = 1+t','y = 2+t','Location','Best')%凡例の設定,'Best'は'SouthWest'や'NorthEast'などに変更可能
    
    saveas(gcf,'output.eps' , 'epsc2'); %epsの色出力
    
end

