/*
離散化＋sontag

sample input
1.00000 1.00000
54
17	34
18	34
19	34
20	34
21	34
22	34
23	34
24	34
25	34
26	34
27	34
28	34
29	34
30	34
31	34
32	34
33	34
34	34
17	17
18	17
19	17
20	17
21	17
22	17
23	17
24	17
25	17
26	17
27	17
28	17
29	17
30	17
31	17
32	17
33	17
34	17
34	17
34	18
34	19
34	20
34	21
34	22
34	23
34	24
34	25
34	26
34	27
34	28
34	29
34	30
34	31
34	32
34	33
34	34
45.00000 45.00000
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#define F "result.csv"
#define Fxy "resultXY.csv"
#define map_size_x 50.0
#define map_size_y 50.0
#define density 4       //離散化の細かさを決める
#define rin 3.0         //全体離散点のノードの半径
#define DT 0.01

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

class WALL{
public:
	double wx, wy;
	WALL(double w_x_, double w_y_) : wx(w_x_), wy(w_y_){}
};

class GRAPH{
public:
	vector<NODE> node;
	vector<STAK> stak;
	vector<vector<EDGE> > edge;
	vector<WALL> wall;
	void print(void);//グラフの情報表示
};

//離散点群を出力．matlabで図を確認．----------------------------------------------
void GRAPH::print(void){
	ofstream outputfile(F);

	for(int i=0; i < node.size(); i++){
		outputfile<< node[i].x<<", "<< node[i].y<<", "<< node[i].r<<endl;
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

//点と点の距離-------------------------------------------------------------------
double point_len(double x, double xx, double y, double yy){
	double min_p_p = 500000.0;
	double length;

	return length = pow( (x - xx) * (x - xx) + (y - yy) * (y - yy), 0.5 );
}

//レイヤー群の最小をとる
double v_min(GRAPH g, double xin, double yin, vector<double> &dist, int hoji_i[1]){
  double v = 0.0;//V(xin, yin)=min_{すべてのノードi}ノードiの関数V(xin,Yin)+始点sからノードiまでの重み
  double min = INF;

	for(int i = 0; i < g.node.size(); i++){
		for(int j = 0; j < g.wall.size(); j++){
			if(point_len(g.node[i].x, g.wall[j].wx, g.node[i].y, g.wall[j].wy) <= g.node[i].r) v = INF;
			else v = g.node[i].v(xin, yin) + dist[i];
		}

		if(v < min){
			min = v;
			hoji_i[0] = i;
		}
	}
	return min;
}

//オイラー法による軌道計算オブジェクト---------------------------------------------
void dissasembled_differential(GRAPH g, double x[2], double p[2], vector<double> &dist){
	int hoji_i[1] = {0};
	int i;
  v_min(g, x[0], x[1], dist, hoji_i);
	i = hoji_i[0];

	double xt = x[0] - g.node[i].x;
	double yt = x[1] - g.node[i].y;
	double sqr = sqrt(xt * xt + yt * yt);

	p[0] = xt / sqr / pow(cos(M_PI / 2.0 / g.node[i].r), 2);
	p[1] = yt / sqr / pow(cos(M_PI / 2.0 / g.node[i].r), 2);
}

// pをuに変換---------------------------------------------------------------------
void euler(GRAPH g, double x[2], double u[2], vector<double> &dist){
	double p[2];
	dissasembled_differential(g, x, p, dist);
	u[0] = -p[0];
	u[1] = -p[1];
}

//オイラー法---------------------------------------------------------------------
void calc_euler(GRAPH g, double start[2], vector<double> &dist, int s){
	double u[2] = {0.0, 0.0};
	double x[2] ={start[0], start[1]};
	double time = 0.0;

	ofstream outputfile_xy(Fxy);
	outputfile_xy << x[0]<< ","<< x[1]<< ","<< 0.0<< endl;
	for (time = 0.0; time <= 200.0 && abs(x[0] - g.node[s].x) > 0.000001 && abs(x[1] - g.node[s].y) > 0.000001; time += DT ){
		euler(g, x, u, dist);
		x[0] += DT * u[0];
		x[1] += DT * u[1];
		cout << x[0]<<" "<<x[1] << endl;
		outputfile_xy << x[0]<< ","<< x[1]<< ","<< time<< endl;
	}
	outputfile_xy << x[0]<< ","<< x[1]<< ","<< time<< endl;
	outputfile_xy.close();
}

//graph描画ファイル生成----------------------------------------------------------
void graph(GRAPH g, vector<double> &dist){
  double xmin = -5.0;
  double xmax = double(map_size_x) + 5.0;
  double ymin = -5.0;
  double ymax = double(map_size_y) + 5.0;
  //double zmin = -1.0;
  double zmax = INF;
  double meshwidth_x = 1.0;
  double meshwidth_y = 1.0;

	FILE *fp;
	double x, y, z;
	int i = 0;
	int dummy[1] = {0};
	fp = fopen("result_3d_z.txt", "w");

	for(y = ymin; y < ymax; y += meshwidth_y){
		for(x = xmin; x < xmax; x += meshwidth_x){
			z = v_min(g, x, y, dist, dummy);
			for(int k = 0; k < g.wall.size(); k++){
				if(g.wall[k].wx - 0.001 <= x && x <= g.wall[k].wx - 0.001
					&& g.wall[k].wy - 0.001 <= y && y <= g.wall[k].wy + 0.001) z = INF;//なくていいかも
			}
			i++;
			cout << i << endl;
			fprintf(fp, "%lf ", z);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("result_3d_x.txt","w");
	for(x = xmin; x < xmax; x += meshwidth_x){
		fprintf(fp, "%lf ", x);
	}
	fclose(fp);
	fp = fopen("result_3d_y.txt","w");
	for(y = ymin; y < ymax; y += meshwidth_y){
		fprintf(fp, "%lf ", y);
	}
	fclose(fp);
}

//edge生成----------------------------------------------------------------------
make_edge(GRAPH *g){
	double node_cost;
	for(int i = 0; i < g->node.size(); i++){
		for(int j = 0; j < g->node.size(); j++){

			double dx = g->node[i].x - g->node[j].x;
      double dy = g->node[i].y - g->node[j].y;
      double nr = g->node[j].r * g->node[j].r;

			if(dx * dx + dy * dy < nr ){
				for(int k = 0; k < g->wall.size(); k++){
					if(point_len(g->node[i].x, g->wall[k].wx, g->node[i].y, g->wall[k].wy) > g->node[i].r
					&& g->node[i].b != -5
					&& g->node[i].r / 2.1 < g->node[j].r
					&& g->node[j].r / 2.1 < g->node[i].r){
						node_cost = double(pow(2, g->node[i].a));
						g->edge[j].push_back(EDGE(i, (g->node[i].v(g->node[j].x, g->node[j].y) + 1000.0 / node_cost)));
					}
					else g->node[i].b = -5;
				}
			}
		}
	}
}
//関数removeの定義---------------------------------------------------------------
template<typename T>//remove(g->node, i)　g->nodeのi番目を削除
void remove(std::vector<T>& vector, unsigned int index){
	vector.erase(vector.begin() + index);
}
//指定した離散点を削除-----------------------------------------------------------
void erase_node(GRAPH *g){
	for(int w = 0; w < g->wall.size(); w++){
		for(int i = 0; i < g->node.size(); i++){
			if(point_len(g->node[i].x, g->wall[w].wx, g->node[i].y, g->wall[w].wy) < 0.001
			 || g->node[i].x < 0.0 || g->node[i].y < 0.0){
				remove(g->node, i);
			}
		}
	}
}
// //もしも壁に近すぎる離散点があれば、ここで削除　　なくていいかも
// void erase_close_node(GRAPH *g){
// 	for(int i = 0; i < g->node.size(); i++){
// 		for(int j = 0; j < g->wall.size(); j++){
// 			if(point_len(g->node[i].x, g->wall[j].wx, g->node[i].y, g->wall[j].wy) < g->node[i].r){
// 				remove(g->node, i);
// 			}
// 		}
// 	}
// }
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
      	&& point_len(dz_x, wall_x, dz_y, wall_y) < s_r
      	&& point_len(dz_x, wall_x, dz_y, wall_y) > s_rr
      	&& r < density - 2){
				if(check_same(g, dz_x, dz_y, s_rr) == 0)
					g->node.push_back(NODE(dz_x, dz_y, s_rr, ain, bin, ccin));
				}
			}
			if(s_rr < g->node[0].r * pow(0.5, density)) break;//十分な細かさが実現したら終了
			else if(point_len(s_xx, wall_x, s_yy, wall_y) < s_rr)//障害物をノード内に含むとき,より細かく離散化．
			make_small_node(g, s_xx, s_yy, s_rr, wall_x, wall_y);
			else{//障害物をノード内に含まないため離散点とする．
				if(check_same(g, s_xx, s_yy, s_rr) == 0)
				g->node.push_back(NODE(s_xx, s_yy, s_rr, ain, bin, ccin));
			}
		}
	}


//平面全体の離散化---------------------------------------------------------------
void desti_all(GRAPH *g){
  double xin, yin;
	int ain = 0, bin = -1, ccin = 0;
  int j = 0;

	for(int x = 0; x * 2 < map_size_x; x++){
  	for(int y = 0; y < map_size_y; y++){
  		if(x % 2 == 0 && y % 4 == 0){
  			xin = y;
				for(int w = 0; w < g->wall.size(); w++){
					if(point_len(xin, g->wall[w].wx, yin, g->wall[w].wy) < rin)
						g->stak.push_back(STAK(xin, yin, rin));
					else{
						if(check_same(g, xin, yin, rin) == 0)
							g->node.push_back(NODE(xin, yin, rin, ain, bin, ccin));
					}
				}
			}
			if(x % 2 == 1 && y % (2 + 4 * j) == 0 && y != 0){
				xin = y;
				for(int w = 0; w < g->wall.size(); w++){
					if(point_len(xin, g->wall[w].wx, yin, g->wall[w].wy) < rin)
						g->stak.push_back(STAK(xin, yin, rin));
					else{
						if(check_same(g, xin, yin, rin) == 0)
							g->node.push_back(NODE(xin, yin, rin, ain, bin, ccin));
					}
				}
				j++;
			}
		}
		yin += 2.0;
  	j = 0;
	}
}

//離散化を管理する関数-----------------------------------------------------------
void discretization(GRAPH *g){
	double s_x, s_y, s_r;
	desti_all(g);
	for(int i = 0; i <= g->stak.size(); i++){
		s_x = g->stak[i].sx;
		s_y = g->stak[i].sy;
		s_r = g->stak[i].sr;
		for(int j=0; j< g->wall.size(); j++){
			if(point_len(s_x, g->wall[j].wx, s_y, g->wall[j].wy));
			make_small_node(g, s_x, s_y, s_r, g->wall[j].wx, g->wall[j].wy);
		}
	}
	erase_node(g);
//	erase_close_node(g);
}

//データインポート--初期地＆壁----------------------------------------------------
void import(GRAPH *g, double start[2]){
  double xin, yin, w_x, w_y;
	int ain = 0, bin = -1, ccin = 0;
	int N;

	cin >> xin >> yin;
	cin >> N;
	start[0] = xin;
	start[1] = yin;

	//壁一個もないときは，(-100, -100)（ロボットが行きえない座標）に壁があるとみなす．
	//g.wall.size()のループがあるため
	if(N == 0)g->wall.push_back(WALL(-100.0, -100.0));

	//壁座標インポート
	for(int n = 0; n < N; n++){
		cin >> w_x >> w_y;
		g->wall.push_back(WALL(w_x, w_y));
	}
}
//点と点の距離でゴールノード番号を取得---------------------------------------------
void goal_import(GRAPH *g, int *s){
	double goal[2];
	double min_p_p = 500000.0;
	double length;

	cin >> goal[0] >> goal[1];

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
  vector<double> dist; // 最短距離

  GRAPH g;
  import(&g, start);//ファイルから入力
	discretization(&g);//空間離散化
	goal_import(&g, &s);//目標座標設定
	g.edge=vector<vector<EDGE> >(g.node.size());
	make_edge(&g);//エッジ作成
	bellman_ford(s, g, dist);//最短経路導出
	//graph(g, dist);//３D画像ファイルに出力
  g.print();//離散点の２D画像ファイルに出力
	calc_euler(g, start, dist, s);//移動経路計算
  return 0;
}
