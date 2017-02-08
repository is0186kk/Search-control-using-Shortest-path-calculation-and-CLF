/*
離散化＋CLF作成
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#define F "result.csv"
#define F3d_x "result_3d_x.txt"
#define F3d_y "result_3d_y.txt"
#define F3d_z "result_3d_z.txt"
#define map_size_xmin 0.0
#define map_size_ymin 0.0
#define map_size_xmax 50.0 //マップサイズ（横）
#define map_size_ymax 50.0 //マップサイズ（縦）
#define detail 4        //離散化の粗さを決める
#define density 4       //離散化の細かさを決める
#define rin 3.0         //全体離散点のノードの半径
#define DT 0.01
#define point_len(x, xx, y, yy)((x - xx) * (x - xx) + (y - yy) * (y - yy))//点と点の距離
using namespace std;

const int INF = 60000;

//class-------------------------------------------------------------------------
class NODE{
	public:
	double x, y, r;
	int a, b, c;
	NODE(double x_, double y_, double r_, int a_, int b_, int c_) : x(x_), y(y_), r(r_), a(a_), b(b_), c(c_){}//( x座標, y座標, 半径, 離散点番号, 壁座標反映後の更新判定に用いる)
	double v(double xin, double yin){
		double xt = xin - x;
		double yt = yin - y;
		double  v = 2.0 * r / M_PI * tan(M_PI / 2.0 / r * sqrt(xt * xt + yt * yt));

		if( xt * xt + yt * yt < r * r && v < INF) return v;
		return INF;
	}

	void print(void){
		cout << x << " " << y << " " << r << endl;
	}
};

class EDGE{
public:
	int to;
	double cost;
	//（どのノードにつながるか, エッジコスト）
	EDGE(int to_, double cost_) : to(to_), cost(cost_) {}
};

//ノード内に障害物を含む場合，この配列に格納．細かな離散化を行う基準点とする．
class STAK{
public:
	double sx, sy, sr;
  STAK(double s_x_, double s_y_, double s_r_) : sx(s_x_), sy(s_y_), sr(s_r_){}
};

class GRAPH{
public:
	vector<NODE> node;
	vector<STAK> stak;
	vector<vector<EDGE> > edge;
	void print(void);//グラフの情報表示
};

//離散点群を出力．matlabで図を確認．----------------------------------------------
void GRAPH::print(void){
	ofstream outputfile(F);

	for(int i=0; i < node.size(); i++){
		outputfile << node[i].x <<", " << node[i].y <<endl;
	}
	outputfile.close();
}

// //グラフの情報表示------------------------------------------------------------
// void GRAPH::print(void){
// 	cout << "number of nodes : " << node.size() << endl;
// 	for(int i = 0; i < node.size(); i++){
// 		node[i].print();
// 	}
// 	cout << "edges : " << endl;
// 	for(int i = 0; i < edge.size(); i++){
// 		for(int t = 0; t < edge[i].size(); t++){
// 			cout << "from " << i << "to " << edge[i][t].to
// 			     << "cost " << edge[i][t].cost << endl;
// 		}
// 	}
// }

//障害物座標データ
class OBSTACLES{
private:
	typedef struct{
		double x;
		double y;
	}xypoint;
public:
	double xmin,xmax,ymin,ymax;
	double rmax;
	vector< vector< vector<xypoint> > > data;//三次元配列に障害物の位置情報を保存
	//特定のグリッドにある点を全て列挙する
	void print_obstacle(int y, int x){
		for(int i = 0; i < data[y][x].size(); i++){
			cout << "obstacle: " << data[y][x][i].x << " " <<data[y][x][i].y << endl;
		}
	}

	//初期化
	OBSTACLES(const char *str,
		double xmin_, double xmax_,
		double ymin_, double ymax_,
		double rmax_){
			xmin = xmin_;
			xmax = xmax_;
			ymin = ymin_;
			ymax = ymax_;
			rmax = rmax_;
			std::cout << "Read obstacle data from "<< str << std::endl;
			std::cout << "x coordinate " << xmin << "--" << xmax << std::endl;
			std::cout << "y coordinate " << ymin << "--" << ymax << std::endl;
			std::cout << "max radius " << rmax << std::endl;
			//行列の要素数を計算
			int xnum = (int)((xmax - xmin) / rmax) + 1;
			int ynum = (int)((ymax - ymin) / rmax) + 1;
			//配列のサイズ変更
			data.resize(ynum);		// ()内の数字が要素数になる
			for(int i = 0; i < ynum; i++ ){
				data[i].resize(xnum);
			}
			//データ読み込み
			std::ifstream ifs(str);
			if (!ifs) {
				cerr << "file failure" << endl;
			}
			while (!ifs.eof()){
				double xin,yin;
				//とりあえず入力
				ifs >> xin >> yin;
				//どの区画に入るか計算
				int y = (int)((yin - ymin) / rmax);
				int x = (int)((xin - xmin) / rmax);
				//当該区画にpush
				if(0 <= y && y < ynum && 0 <= x && x < xnum){
					xypoint inpoint={xin, yin};
					data[y][x].push_back(inpoint);
				}
			}
	}
};

//bellman_ford法で経路コスト計算--------------------------------------------------
bool bellman_ford(int s, GRAPH g, vector<double> &dist){//nは頂点数、sは開始頂点
  dist = vector<double>(g.node.size(), INF);
  dist[s] = 0; // 開始点の距離は0
  for (int i = 0; i < g.node.size(); i++) {
    for (int v = 0; v < g.node.size(); v++) {
      for (int k = 0; k < g.edge[v].size(); k++) {
			  EDGE e = g.edge[v][k];
        if (dist[v] != INF && dist[e.to] > dist[v] + e.cost){
          dist[e.to] = dist[v] + e.cost;
          if (i == g.node.size() - 1) return true; //n回目にも更新があるなら負の閉路が存在
        }
      }
	  }
  }
  return false;
}
// 配列内に同じ座標が存在する場合にpush_backしないための関数------------------------
// 同じ座標に離散点を作ろうとしていた場合は半径が大きいほうを削除
int check_same(GRAPH *g, double x, double y, double r){
	for(int i = 0; i < g->node.size(); i++){
		if(g->node[i].x - 0.001 <= x && x <= g->node[i].x + 0.001
			&& g->node[i].y - 0.001 <= y && y <= g->node[i].y + 0.001){
			if(g->node[i].r > r){
				g->node[i].r = r;
				return 0;
			}
			else return 1;
		}
	}
	return 0;
}

//障害物近傍の大きいノードを離散化-------------------------------------------------
void make_small_node(GRAPH *g, double s_x, double s_y, double s_r, double wall_x, double wall_y){
	double s_xx, s_yy, s_rr;
	int ain = 0, bin = -1, ccin = 0;
  double r = log2(g->node[0].r / s_r);
  const double dx[5] = {0, -1, +1, -1, +1};
  const double dy[5] = {0, +1, +1, -1, -1};
  const double dz[4][2] = {{+1, 0}, {-1, 0}, {0, +1}, {0, -1}};
  double dz_x, dz_y;
	ain = r;

	//中心，左上，右上，左下，右下　の順に大きい円の上に小さい円を並べる．
	for(int i = 0; i < 5; i++){
		s_rr = s_r * 0.5;
    s_xx = s_x + dx[i] * pow(0.5, r);
    s_yy = s_y + dy[i] * pow(0.5, r);

  //障害物を境に左上点群,右上点群,左下点群,右下点群ができるので,群同士を繋ぐ離散点を設計　
  //離散幅が細かすぎず,障害物にも近すぎないとき,中央座標から一定距離の点を離散点とする．
		for(int j = 0; j < 4; j++){
    	dz_x = s_x + dz[j][0] * 2.0 * pow(0.5, r);
    	dz_y = s_y + dz[j][1] * 2.0 * pow(0.5, r);
			if(i == 0
      	&& point_len(dz_x, wall_x, dz_y, wall_y) < s_r * s_r
      	&& point_len(dz_x, wall_x, dz_y, wall_y) > s_rr * s_rr
      	&& r < density - 1){
				if(check_same(g, dz_x, dz_y, s_rr) == 0)
					g->node.push_back(NODE(dz_x, dz_y, s_rr, ain, bin, ccin));
			}
		}
		if(s_rr < g->node[0].r * pow(0.5, density)) break;//十分な細かさが実現したら終了
		else if(point_len(s_xx, wall_x, s_yy, wall_y) < s_rr * s_rr)//障害物をノード内に含むとき,より細かく離散化．
		make_small_node(g, s_xx, s_yy, s_rr, wall_x, wall_y);
		else{//障害物をノード内に含まないため離散点とする．
			if(check_same(g, s_xx, s_yy, s_rr) == 0)
			g->node.push_back(NODE(s_xx, s_yy, s_rr, ain, bin, ccin));
		}
	}
}

//壁と近いとき1を返す
int check_wall(GRAPH* g, OBSTACLES obstacle, double xin, double  yin, double rrin, int a){
	int ynum = (int)((yin - obstacle.ymin) / obstacle.rmax);
	int xnum = (int)((xin - obstacle.xmin) / obstacle.rmax);

	int ynum_min = ynum - 1;
	int ynum_max = ynum + 1;
	int xnum_min = xnum - 1;
	int xnum_max = xnum + 1;

	if(ynum_min < 0)ynum_min = 0;
	if(xnum_min < 0)xnum_min = 0;
	if(ynum_max > (int)((obstacle.ymax - obstacle.ymin) / obstacle.rmax))
		 ynum_max = (int)((obstacle.ymax - obstacle.ymin) / obstacle.rmax);
	if(xnum_max > (int)((obstacle.xmax - obstacle.xmin) / obstacle.rmax))
		 xnum_max = (int)((obstacle.xmax - obstacle.xmin) / obstacle.rmax);

	for(int j = ynum_min; j < ynum_max; j++){
		for(int i = xnum_min; i < xnum_max; i++){
			for(int k = 0; k < obstacle.data[j][i].size(); k++){
				if(point_len(xin, obstacle.data[j][i][k].x, yin, obstacle.data[j][i][k].y) < rrin * rrin){
					if(a == 0)return 1;
					else make_small_node(g, xin, yin, rrin, obstacle.data[j][i][k].x, obstacle.data[j][i][k].y);
				}
			}
		}
	}
	return 0;
}

//レイヤー群の最小をとる----------------------------------------------------------
double v_min(GRAPH *g, OBSTACLES obstacle,
						 double xin, double yin, vector<double> &dist, int hoji_i[1]){
  double v = 0.0;//V(xin, yin)=min_{すべてのノードi}ノードiの関数V(xin,Yin)+始点sからノードiまでの重み
  double min = INF;
	int ynum = (int)((yin - obstacle.ymin) / obstacle.rmax);
	int xnum = (int)((xin - obstacle.xmin) / obstacle.rmax);

	int ynum_min = ynum - 1;
	int ynum_max = ynum + 1;
	int xnum_min = xnum - 1;
	int xnum_max = xnum + 1;

	if(ynum_min < 0)ynum_min = 0;
	if(xnum_min < 0)xnum_min = 0;
	if(ynum_max > (int)((obstacle.ymax - obstacle.ymin) / obstacle.rmax))
		 ynum_max = (int)((obstacle.ymax - obstacle.ymin) / obstacle.rmax);
	if(xnum_max > (int)((obstacle.xmax - obstacle.xmin) / obstacle.rmax))
		 xnum_max = (int)((obstacle.xmax - obstacle.xmin) / obstacle.rmax);

	for(int p = 0; p < g->node.size(); p++){
		v = g->node[p].v(xin, yin) + dist[p];
		for(int j = ynum_min; j < ynum_max; j++){
			for(int i = xnum_min; i < xnum_max; i++){
				for(int k = 0; k < obstacle.data[j][i].size(); k++){
					if(point_len(g->node[p].x, obstacle.data[j][i][k].x, g->node[p].y, obstacle.data[j][i][k].y) <= g->node[p].r * g->node[p].r)
					v = INF;
					//else v = g->node[p].v(xin, yin) + dist[p];
					}
				}
			}
		if(v < min){
			min = v;
			hoji_i[0] = p;
		}
	}
	return min;
}

//graph描画ファイル生成----------------------------------------------------------
void graph(GRAPH *g, OBSTACLES obstacle, vector<double> &dist){
  double xmin = double(map_size_xmin) - 5.0;
  double xmax = double(map_size_xmax) + 5.0;
  double ymin = double(map_size_ymin) - 5.0;
  double ymax = double(map_size_ymax) + 5.0;
  double meshwidth_x = 1.0;
  double meshwidth_y = 1.0;

	double x, y, z;
	int i = 0;
	int dummy[1] = {0};

	ofstream outputfile_z(F3d_z);
	for(y = ymin; y <= ymax; y += meshwidth_y){
		for(x = xmin; x <= xmax; x += meshwidth_x){
			z = v_min(g, obstacle, x, y, dist, dummy);
			int ynum = (int)((y - obstacle.ymin) / obstacle.rmax);
			int xnum = (int)((x - obstacle.xmin) / obstacle.rmax);
			for(int k = 0; k < obstacle.data[ynum][xnum].size(); k++){
			 	if(point_len(x, obstacle.data[ynum][xnum][k].x,
			 		           y, obstacle.data[ynum][xnum][k].y) < 0.01) z = INF;
			}
			i++;
			cout << i << endl;
			outputfile_z << z << " ";
		}
		outputfile_z << endl;
	}
	outputfile_z.close();

	ofstream outputfile_x(F3d_x);
	for(x = xmin; x <= xmax; x += meshwidth_x){
		outputfile_x << x << " ";
	}
	outputfile_x.close();
	ofstream outputfile_y(F3d_y);
	for(y = ymin; y <= ymax; y += meshwidth_y){
		outputfile_y << y << " ";
	}
	outputfile_y.close();
}

//関数removeの定義---------------------------------------------------------------
template<typename T>//remove(g->node, i)　g->nodeのi番目を削除
void remove(std::vector<T>& vector, unsigned int index){
	vector.erase(vector.begin() + index);
}
//指定した離散点を削除------------------------------------------------------------
//壁に近すぎる離散点を消す
void erase_wall(GRAPH *g, OBSTACLES obstacle){
	for(int i = 0; i < g->node.size(); i++){
		int y = (int)((g->node[i].y - obstacle.ymin) / obstacle.rmax);
		int x = (int)((g->node[i].x - obstacle.xmin) / obstacle.rmax);

		for(int w = 0; w < obstacle.data[y][x].size(); w++){
			if(point_len(g->node[i].x, obstacle.data[y][x][w].x,
				           g->node[i].y, obstacle.data[y][x][w].y) < 0.001){
									 remove(g->node, i);
									 i--;
			}
		}
	}
}

//平面全体の離散化
void desti_all2(GRAPH *g, OBSTACLES obstacle, double xin, double yin){
	int ain = 0, bin = -1, ccin = 0;
	if(check_wall(g, obstacle, xin, yin, rin, 0) == 0){
		if(check_same(g, xin, yin, rin) == 0)
		g->node.push_back(NODE(xin, yin, rin, ain, bin, ccin));
	}
	else g->stak.push_back(STAK(xin, yin, rin));
}

//平面全体の離散化---------------------------------------------------------------
void desti_all(GRAPH *g, OBSTACLES obstacle){
	int size_x = map_size_xmax / detail + 1.0;
	int size_y = map_size_ymax / detail + 1.0;

	for(int y = 1; y < size_y; y++){
		for(int x = 1; x < size_x; x++){
			double xin = detail * x;
			double yin = detail * y;
			desti_all2(g, obstacle, xin, yin);
			if((x - 1.0 / 2.0) < size_x  || (y - 1.0 / 2.0) < size_y){
				xin = detail * (x - 1.0 / 2.0);
				yin = detail * (y - 1.0 / 2.0);
			}
			desti_all2(g, obstacle, xin, yin);
		}
	}
}

//離散化を管理する関数-----------------------------------------------------------
void discretization(GRAPH *g, OBSTACLES obstacle){
	desti_all(g, obstacle);
	for(int s = 0; s < g->stak.size(); s++){
		check_wall(g, obstacle, g->stak[s].sx, g->stak[s].sy, g->stak[s].sr, 1);
	}
	cout << g->node.size() << endl;
	erase_wall(g, obstacle);
	cout << g->node.size() << endl;
}

//edge生成--------------------------------------------------------------------
make_edge(GRAPH *g, OBSTACLES obstacle){
	for(int i = 0; i < g->node.size(); i++){
		for(int j = 0; j < g->node.size(); j++){
			double dx = g->node[i].x - g->node[j].x;
      double dy = g->node[i].y - g->node[j].y;
      double nr = g->node[j].r * g->node[j].r;

			if(dx * dx + dy * dy < nr){
				if(check_wall(g, obstacle, g->node[i].x, g->node[i].y, g->node[i].r, 0) == 0
				&& g->node[i].b != -5
				&& g->node[i].r / 2.1 < g->node[j].r
				&& g->node[j].r / 2.1 < g->node[i].r){
					int node_cost = pow(2, g->node[i].a);
					g->edge[j].push_back(EDGE(i, (g->node[i].v(g->node[j].x, g->node[j].y) + 1000.0 / double(node_cost))));
				}
				else g->node[i].b = -5;
			}
		}
	}
}

//データインポート--初期地＆壁----------------------------------------------------
void import(GRAPH *g, double start[2]){
  double xin, yin, w_x, w_y;
	int N;
	xin = 5.0;
	yin = 5.0;
	start[0] = xin;
	start[1] = yin;
}
//点と点の距離でゴールノード番号を取得---------------------------------------------
void goal_import(GRAPH *g, int *s){
	double goal[2];
	double min_p_p = INF;
	double length;
	double xin, yin;
	xin = 49;
	yin = 49;

	goal[0]=xin;
	goal[1]=yin;

	for(int i = 0; i < g->node.size(); i++){
		length = point_len(goal[0], g->node[i].x, goal[1], g->node[i].y);
		if(length < min_p_p){
			*s = i;
			min_p_p = length;
		}
	 }
}
//main--------------------------------------------------------------------------
int main() {
  int s;
	double start[2];
	double xmin = 0;
	double xmax = map_size_xmax;
	double ymin = 0;
	double ymax = map_size_ymax;
	double rmax = 10;

  vector<double> dist; // 最短距離
  GRAPH g;
  import(&g, start);//ファイルから入力
	OBSTACLES obstacle("test_dots.txt", xmin, xmax, ymin, ymax, rmax);
	//obstacle.print_obstacle(5, 3);
	discretization(&g, obstacle);//空間離散化
	goal_import(&g, &s);//目標座標設定
	g.edge=vector<vector<EDGE> >(g.node.size());
	make_edge(&g, obstacle);//エッジ作成
	bellman_ford(s, g, dist);//最短経路導出
	graph(&g, obstacle, dist);//３D画像ファイルに出力
  g.print();//離散点の２D画像ファイルに出力
  return 0;
}
