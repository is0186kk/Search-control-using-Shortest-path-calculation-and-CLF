/*
離散化＋sontag

入力フォーマット

＝＝＝シミュレーション設定を記述するテキストファイルの内容＝＝＝
障害物情報を記述するテキストファイルのファイル名
ロボットの初期位置のｘ座標　ｙ座標
目標位置のｘ座標　ｙ座標

＝＝＝シミュレーション設定を記述するテキストファイルの内容＝＝＝

＝＝＝障害物情報を記述するテキストファイルの内容＝＝＝
map_size_xmin map_size_ymin
map_size_xmax map_size_ymax
障害物1のｘ座標　障害物1のｙ座標
障害物2のｘ座標　障害物2のｙ座標
・
・
・
＝＝＝障害物情報を記述するテキストファイルの内容ここまで＝＝＝


*/
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <float.h>//double型の最大値が保存されている

using namespace std;
#define point_len(x, xx, y, yy)((x - xx) * (x - xx) + (y - yy) * (y - yy))//点と点の距離
const double INF = DBL_MAX;             //double型の最大値

//------------------------------------------------------------------------------------------------------
//シミュレーション条件に関する設定
//------------------------------------------------------------------------------------------------------
//const char *simulation_condition = "condition_lab.txt";//シミュレーション条件を記述してあるファイル名
const char *simulation_condition = "condition_simple.txt";//シミュレーション条件を記述してあるファイル名

//------------------------------------------------------------------------------------------------------
//空間の離散化に関する設定
//------------------------------------------------------------------------------------------------------
char obstacle_file[256];               //障害物データが記録されているテキストファイルのファイル名
const char *Fname_matlab = "result.csv";//離散点群をmatlabのフォーマットで出力する際のファイル名
const double detail = 4;                //離散化の第1レベルにおける粗さ(点をどのくらいの幅で配置するか)
const double rin = 3.0;                 //第1レベルの離散化における円の半径.rin > sqrt(2.0)/2.0 * detailを満たす必要あり
const int density = 2;                  //離散化の再起プロセスを何回くりかえすか
const double K = (int)(rin + 0.5);        //K>rinである必要がある
//------------------------------------------------------------------------------------------------------
//オイラー法に関するパラメータ
//------------------------------------------------------------------------------------------------------
const char *Fxy = "resultXY.csv";       //ログの出力先
const double euler_end_time = 100.0;    //シミュレーションの最大時間
const double DT = 0.01;                 //オイラー法の時間パラメータ
double start[2];                        //初期状態量
//------------------------------------------------------------------------------------------------------
//制御目標に関するパラメータ
//------------------------------------------------------------------------------------------------------
double goal[2];                         //目標位置

//class-------------------------------------------------------------------------
//a : 離散化を細かくするかの敷居値，０は全体の離散化（レベル１）　１はレベル２以上
//b :
//c :
class NODE{
public:
	double x, y, r;
	int a, b, c;
	NODE(double x_, double y_, double r_, int a_, int b_) : x(x_), y(y_), r(r_), a(a_), b(b_){}//( x座標, y座標, 半径, 離散点番号, 壁座標反映後の更新判定に用いる)
	double v(double xin, double yin){
		double xt = xin - x;
		double yt = yin - y;
		double rt = sqrt(xt * xt + yt * yt);
		if (rt < r){
			double  v = 2.0 * r / M_PI * tan(M_PI / 2.0 / r * rt);
			return v;
		}
		else
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


class GRAPH{
public:
	vector<NODE> node;
	vector<vector<EDGE> > edge;
	//与えられた座標に一番近い中心を持つノード番号を返す
	void calc_nodenumber(const double x, const double y, int *s){
		double min_p_p = INF;
		double length;

		//cin >> goal[0] >> goal[1];
		for (int i = 0; i < node.size(); i++){
			length = point_len(x, node[i].x, y, node[i].y);
			if (length < min_p_p){
				*s = i;
				min_p_p = length;
			}
		}
	}

	void print(void);//グラフの情報表示
};

//離散点群を出力．matlabで図を確認．----------------------------------------------
void GRAPH::print(void){
	ofstream outputfile(Fname_matlab);

	for (int i = 0; i < node.size(); i++){
		outputfile << node[i].x << ", " << node[i].y << ", " << node[i].r << ", " << node[i].a << endl;
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


class OBSTACLES{
private:
	typedef struct{
		double x;
		double y;
	}xypoint;
public:
	double xmin, xmax, ymin, ymax;
	double lwidth;
	vector< vector< vector<xypoint> > > data;//三次元配列に障害物の位置情報を保存

	//特定のグリッドにある点を全て列挙する
	void print_obstacle(int y, int x){
		for (int i = 0; i < data[y][x].size(); i++){
			cout << "obstacle: " << data[y][x][i].x << " " << data[y][x][i].y << endl;
		}
	}

	//初期化
	OBSTACLES(const char *str){
		lwidth = K;
		//データ読み込み準備
		std::ifstream ifs(str);
		if (!ifs) {
			cerr << "file failure" << endl;
		}
		//ヘッダ読み込み
		ifs >> xmin >> xmax;
		ifs >> ymin >> ymax;

		//データ本体読み込み準備
		std::cout << "Read obstacle data from " << str << std::endl;
		std::cout << "x coordinate " << xmin << "--" << xmax << std::endl;
		std::cout << "y coordinate " << ymin << "--" << ymax << std::endl;
		std::cout << "koushi haba " << lwidth << std::endl;
		//行列の要素数を計算
		int xnum = (int)((xmax - xmin) / lwidth) + 1;
		int ynum = (int)((ymax - ymin) / lwidth) + 1;
		//配列のサイズ変更
		data.resize(ynum);		// ()内の数字が要素数になる
		for (int i = 0; i < ynum; i++){
			data[i].resize(xnum);
		}

		//データ本体読み込み
		while (!ifs.eof()){
			double xin, yin;
			//とりあえず入力
			ifs >> xin >> yin;
			//どの区画に入るか計算
			int y = (int)((yin - ymin) / lwidth);
			int x = (int)((xin - xmin) / lwidth);
			//当該区画にpush
			if (0 <= y && y < ynum && 0 <= x && x < xnum){
				xypoint inpoint = { xin, yin };
				data[y][x].push_back(inpoint);
			}
		}
	}
	//与えられた（ｘ，ｙ）座標が属する区画の九近傍（前後左右斜め）の区画座標を返す関数
	void calc_neighborhood(double xin, double yin, int *xoutmin, int *xoutmax, int *youtmin, int *youtmax){
		int ynum = (int)((yin - ymin) / lwidth);
		int xnum = (int)((xin - xmin) / lwidth);

		int ynum_min = ynum - 1;
		int ynum_max = ynum + 1;
		int xnum_min = xnum - 1;
		int xnum_max = xnum + 1;

		double ynum_max_limit = (int)((ymax - ymin) / lwidth);
		double xnum_max_limit = (int)((xmax - xmin) / lwidth);

		if (ynum_min < 0) ynum_min = 0;
		if (xnum_min < 0) xnum_min = 0;
		if (ynum_max > ynum_max_limit) ynum_max = ynum_max_limit;
		if (xnum_max > xnum_max_limit) xnum_max = xnum_max_limit;

		*xoutmin = xnum_min;
		*xoutmax = xnum_max;
		*youtmin = ynum_min;
		*youtmax = ynum_max;
	}
};


//bellman_ford法で経路コスト計算--------------------------------------------------
bool bellman_ford(int s, GRAPH g, vector<double> &dist){//nは頂点数、sは開始頂点
	int g_node_size = g.node.size();
	dist = vector<double>(g_node_size, INF);
	dist[s] = 0; // 開始点の距離は0
	for (int i = 0; i < g_node_size; i++) {
		for (int v = 0; v < g_node_size; v++) {
			int g_edgev_size = g.edge[v].size();
			for (int k = 0; k < g_edgev_size; k++) {
				EDGE e = g.edge[v][k];
				if (dist[v] != INF && dist[e.to] > dist[v] + e.cost){
					dist[e.to] = dist[v] + e.cost;
					if (i == g_node_size - 1) return true; //n回目にも更新があるなら負の閉路が存在
				}
			}
		}
	}
	return false;
}

//壁と近いとき1を返す．また，obstaclex, obstacleyに近い壁の座標を出力する．
//出力される壁は一番最初に見つかったもので，最近近傍のものとは限らない．
int check_wall(GRAPH* g, OBSTACLES obstacle, double xin, double  yin, double rrin){
	int xnum_min, xnum_max;
	int ynum_min, ynum_max;
	obstacle.calc_neighborhood(xin, yin, &xnum_min, &xnum_max, &ynum_min, &ynum_max);//9近傍を列挙

	for (int j = ynum_min; j < ynum_max; j++){
		for (int i = xnum_min; i < xnum_max; i++){
			for (int k = 0; k < obstacle.data[j][i].size(); k++){
				double obx = obstacle.data[j][i][k].x;
				double oby = obstacle.data[j][i][k].y;
				if (point_len(xin, obx, yin, oby) < rrin * rrin){
					return 1;
				}
			}
		}
	}
	return 0;
}

//レイヤー群の最小をとる----------------------------------------------------------
double v_min(GRAPH *g, OBSTACLES obstacle,
	double xin, double yin, vector<double> &dist, int *hoji_i){
	double v;//V(xin, yin)=min_{すべてのノードi}ノードiの関数V(xin,Yin)+始点sからノードiまでの重み
	double min = INF;

	for (int p = 0; p < g->node.size(); p++){
		double ox = g->node[p].x;
		double oy = g->node[p].y;
		double obr = g->node[p].r;
		if (point_len(xin, ox, yin, oy) >= obr * obr){
			v = INF;
		}
		else{
			v = g->node[p].v(xin, yin) + dist[p];
		}
		if (v < min){
			min = v;
			*hoji_i = p;
		}
	}
	return min;
}
//オイラー法による軌道計算オブジェクト---------------------------------------------
void dissasembled_differential(GRAPH *g, OBSTACLES obstacle, double x[2], double p[2], vector<double> &dist){
	int i;
	v_min(g, obstacle, x[0], x[1], dist, &i);

	double xt = x[0] - g->node[i].x;
	double yt = x[1] - g->node[i].y;
	double sqr = sqrt(xt * xt + yt * yt);
	double k = 1.0 / sqr / pow(cos(M_PI / 2.0 / g->node[i].r), 2);

	p[0] = xt * k;
	p[1] = yt * k;
}

// pをuに変換---------------------------------------------------------------------
void euler(GRAPH *g, OBSTACLES obstacle, double x[2], double u[2], vector<double> &dist){
	double p[2];
	dissasembled_differential(g, obstacle, x, p, dist);
	u[0] = -p[0];
	u[1] = -p[1];
}

//オイラー法本体---------------------------------------------------------------------
void calc_euler(GRAPH *g, OBSTACLES obstacle, vector<double> &dist, int s){
	double u[2] = { 0.0, 0.0 };
	double x[2] = { start[0], start[1] };

	ofstream outputfile_xy(Fxy);

	for (double time = 0.0; time <= euler_end_time; time += DT){
		euler(g, obstacle, x, u, dist);

		cout << x[0] << " " << x[1] << endl;
		outputfile_xy << x[0] << "," << x[1] << "," << time << endl;

		x[0] += DT * u[0];
		x[1] += DT * u[1];
	}

	outputfile_xy.close();
}




//edge生成--------------------------------------------------------------------
void make_edge(GRAPH *g, OBSTACLES obstacle){
	int g_node_size = g->node.size();
	g->edge = vector<vector<EDGE> >(g_node_size);
	for (int i = 0; i < g_node_size; i++){//From i To jにエッジがはれるか
		for (int j = 0; j < g_node_size; j++){
			double dx = g->node[i].x - g->node[j].x;
			double dy = g->node[i].y - g->node[j].y;
			double nr = g->node[j].r * g->node[j].r;

			if (dx * dx + dy * dy < nr){//ノードjの円にノードiの中心があるなら
				int node_cost = pow(2, g->node[j].a);
				double v = g->node[j].v(g->node[i].x, g->node[i].y);
				double edge_cost = (v + 1000.0 / node_cost);
				g->edge[i].push_back(EDGE(j, edge_cost));
			}
		}
	}
}

//条件ファイルよりシミュレーション条件を読み込み------------------------------------
void read_simulation_condition(void){
	std::ifstream ifs(simulation_condition);
	cout << simulation_condition << endl;
	if (!ifs) {
		cerr << "simulation file failure" << endl;
	}
	ifs.getline(obstacle_file, 256 - 1);//障害物データ
	ifs >> start[0] >> start[1];//初期位置
	ifs >> goal[0] >> goal[1];//目標位置
	cout << obstacle_file << endl;
}
void make_node_sub(GRAPH *g, OBSTACLES obstacle, double xin, double yin, double r, int level){
	int iswall = check_wall(g, obstacle, xin, yin, r);
	if (iswall == 0){
		//障害物とかぶってなければノード追加
		g->node.push_back(NODE(xin, yin, r, level, 0));
	}
	else{
		if (level < density){
			double r_next = r / sqrt(2.0) * 1.1;//1.1はマージン
			double d = r / 2.0;
			make_node_sub(g, obstacle, xin + d, yin + d, r_next, level + 1);
			make_node_sub(g, obstacle, xin - d, yin + d, r_next, level + 1);
			make_node_sub(g, obstacle, xin + d, yin - d, r_next, level + 1);
			make_node_sub(g, obstacle, xin - d, yin - d, r_next, level + 1);
		}
	}
}
void make_node(GRAPH *g, OBSTACLES obstacle){
	double d = detail;//この幅で代表点を配置
	double r = d * 3.0 / 4.0;//最初のレベルの円の半径
	int size_x = obstacle.xmax / d + 1.0;
	int size_y = obstacle.ymax / d + 1.0;

	for (int y = 1; y < size_y; y++){
		for (int x = 1; x < size_x; x++){
			double xin = detail * x;
			double yin = detail * y;
			//xin,yin,rの円を配置
			make_node_sub(g, obstacle, xin, yin, r, 0);
			if ((x - 1.0 / 2.0) < size_x || (y - 1.0 / 2.0) < size_y){
				xin = detail * (x - 1.0 / 2.0);
				yin = detail * (y - 1.0 / 2.0);
			}
			make_node_sub(g, obstacle, xin, yin, r, 0);
		}
	}
}
//main--------------------------------------------------------------------------
int main() {
	read_simulation_condition();        //条件ファイルよりシミュレーション条件読み込み

	int s;//ゴールノード番号
	vector<double> dist; // 最短距離
	OBSTACLES obstacle(obstacle_file);
	GRAPH g;

	cout << "makenode start" << endl;
	make_node(&g, obstacle);
	//	discretization(&g, obstacle);       //空間離散化によるノード作成
	cout << "nodenum : " << g.node.size() << endl;
	cout << "makenode end" << endl;

	cout << "makeedge start" << endl;
	make_edge(&g, obstacle);            //エッジ作成
	cout << "makeedge end" << endl;

	g.calc_nodenumber(goal[0], goal[1], &s);  //与えられた座標に一番近い中心を持つノード番号を計算

	cout << "bellmanford start" << endl;
	bellman_ford(s, g, dist);                 //最短経路計算
	cout << "bellmanford end" << endl;

	g.print();                                //離散点の２D画像ファイルに出力
	calc_euler(&g, obstacle, dist, s);        //移動経路計算
	return 0;
}
