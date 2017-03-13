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
#include <queue>//Dijkstraアルゴリズムのために必要

using namespace std;
#define point_len(x, xx, y, yy)((x - xx) * (x - xx) + (y - yy) * (y - yy))//点と点の距離
const double INF = DBL_MAX;             //double型の最大値

//------------------------------------------------------------------------------------------------------
//シミュレーション条件に関する設定
//------------------------------------------------------------------------------------------------------
const char *simulation_condition = "condition_lab.txt";//シミュレーション条件を記述してあるファイル名
//const char *simulation_condition = "condition_simple.txt";//シミュレーション条件を記述してあるファイル名

//------------------------------------------------------------------------------------------------------
//オイラー法に関するパラメータ
//------------------------------------------------------------------------------------------------------
const char *Fxy = "resultXY.csv";       //ログの出力先
const double euler_end_time = 400.0;    //シミュレーションの最大時間
const double DT = 0.0001;               //オイラー法の時間パラメータ
double start[2];                        //初期状態量

//------------------------------------------------------------------------------------------------------
//制御目標に関するパラメータ
//------------------------------------------------------------------------------------------------------
double goal[2];                         //目標位置

//------------------------------------------------------------------------------------------------------
//空間の離散化に関する設定
//------------------------------------------------------------------------------------------------------
char obstacle_file[256];                //障害物データが記録されているテキストファイルのファイル名
const char *Fname_matlab = "result.csv";//離散点群をmatlabのフォーマットで出力する際のファイル名
const double discretization_width = 16;  //離散化の第1レベルにおける粗さ(点をどのくらいの幅で配置するか)
const double r_first_level = 12.0;       //第1レベルの離散化における円の半径.
//r_first_level > sqrt(2.0)/2.0 * discretization_widthを満たす必要あり
const int level_max = 3;                //離散化の再起プロセスを何回くりかえすか

//------------------------------------------------------------------------------------------------------
//条件ファイルよりシミュレーション条件を読み込み
//------------------------------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------------------------------
//障害物を取り扱うためのクラス
//------------------------------------------------------------------------------------------------------
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
	//初期化
	OBSTACLES(const char *str);
	//与えられた（ｘ，ｙ）座標が属する区画の九近傍（前後左右斜め）の区画座標を返す関数
	//const修飾子を追加していて，ざっくり言えば，「内部変数をいじるな！」と指定してある．
	void calc_neighborhood(double xin, double yin, int *xoutmin, int *xoutmax, int *youtmin, int *youtmax) const;
	//特定のグリッドにある点を全て列挙する
	void print_obstacle(int y, int x) const;

};

OBSTACLES::OBSTACLES(const char *str){
	lwidth = r_first_level;
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
//const修飾子を追加していて，ざっくり言えば，「内部変数をいじるな！」と指定してある．
void OBSTACLES::calc_neighborhood(double xin, double yin, int *xoutmin, int *xoutmax, int *youtmin, int *youtmax) const{
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
//特定のグリッドにある点を全て列挙する
void OBSTACLES::print_obstacle(int y, int x) const{
	for (int i = 0; i < data[y][x].size(); i++){
		cout << "obstacle: " << data[y][x][i].x << " " << data[y][x][i].y << endl;
	}
}
//------------------------------------------------------------------------------------------------------
//提案法におけるグラフ構造を取り扱うためのクラス
//------------------------------------------------------------------------------------------------------
//ノードの定義-------------------------------------------------------------------
class NODE{
public:
	double x, y, r;//( x座標, y座標, 半径,
	int level;//何レベル目の再起アルゴリズムで追加されたか

	NODE(double x_, double y_, double r_, int level_) : x(x_), y(y_), r(r_), level(level_){}
	double v(double xin, double yin);
	void print(void){
		cout << x << " " << y << " " << r << endl;
	}
};

double NODE::v(double xin, double yin){
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

//エッジの定義-------------------------------------------------------------------
class EDGE{
public:
	int to;
	double cost;
	//（どのノードにつながるか, エッジコスト）
	EDGE(int to_, double cost_) : to(to_), cost(cost_) {}
};

//グラフの定義-------------------------------------------------------------------
class GRAPH{
public:
	//グラフデータ構造
	vector<NODE> node;
	vector<vector<EDGE> > edge;
	//クイックアクセスのためのデータ構造
	double xmin, xmax, ymin, ymax;
	double lwidth;
	int xnum, ynum;
	vector< vector< vector<int> > > data;//三次元配列にノードのインデックスを保存

	//初期化
	GRAPH(OBSTACLES &obstacle);
	//ノード追加
	void add_node(NODE n);
	//与えられた（ｘ，ｙ）座標が属する区画の九近傍（前後左右斜め）の区画座標を返す関数
	//const修飾子を追加していて，ざっくり言えば，「内部変数をいじるな！」と指定してある．
	void calc_neighborhood(double xin, double yin, int *xoutmin, int *xoutmax, int *youtmin, int *youtmax) const;
	//与えられた座標に一番近い中心を持つノード番号を返す
	void calc_nodenumber(const double x, const double y, int *s);
	void print(void);//グラフの情報表示
	void print_full(void);//グラフの情報表示
};
//初期化------------------------------------------------------------------------
GRAPH::GRAPH(OBSTACLES &obstacle){
	xmin = obstacle.xmin;
	xmax = obstacle.xmax;
	ymin = obstacle.ymin;
	ymax = obstacle.ymax;
	lwidth = obstacle.lwidth;
	//行列の要素数を計算
	xnum = (int)((xmax - xmin) / lwidth) + 1;
	ynum = (int)((ymax - ymin) / lwidth) + 1;
	//配列のサイズ変更
	data.resize(ynum);		// ()内の数字が要素数になる
	for (int i = 0; i < ynum; i++){
		data[i].resize(xnum);
	}
}
//ノード追加--------------------------------------------------------------------
void GRAPH::add_node(NODE n){
	int current_index = node.size();//追加するものが何番目のindexか
	node.push_back(n);//データにＰＵＳＨ
	//どの区画に入るか計算
	int y = (int)((n.y - ymin) / lwidth);
	int x = (int)((n.x - xmin) / lwidth);
	//当該区画にノードのindexをpush
	if (0 <= y && y < ynum && 0 <= x && x < xnum){
		data[y][x].push_back(current_index);
	}
}
//与えられた（ｘ，ｙ）座標が属する区画の九近傍（前後左右斜め）の区画座標を返す関数
//const修飾子を追加していて，ざっくり言えば，「内部変数をいじるな！」と指定してある．
void GRAPH::calc_neighborhood(double xin, double yin, int *xoutmin, int *xoutmax, int *youtmin, int *youtmax) const{
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
//与えられた座標に一番近い中心を持つノード番号を返す--------------------------------
void GRAPH::calc_nodenumber(const double x, const double y, int *s){
	double min_p_p = INF;
	double length;

	for (int i = 0; i < node.size(); i++){
		length = point_len(x, node[i].x, y, node[i].y);
		if (length < min_p_p){
			*s = i;
			min_p_p = length;
		}
	}
}
//離散点群を出力．matlabで図を確認．----------------------------------------------
void GRAPH::print(void){
	ofstream outputfile(Fname_matlab);

	for (int i = 0; i < node.size(); i++){
		outputfile << node[i].x << ", " << node[i].y << ", " << node[i].r << ", " << node[i].level << endl;
	}
	outputfile.close();
}
//グラフの情報表示------------------------------------------------------------
void GRAPH::print_full(void){
	cout << "number of nodes : " << node.size() << endl;
	for (int i = 0; i < node.size(); i++){
		node[i].print();
	}
	cout << "edges : " << endl;
	for (int i = 0; i < edge.size(); i++){
		for (int t = 0; t < edge[i].size(); t++){
			cout << "from " << i << "to " << edge[i][t].to
				<< "cost " << edge[i][t].cost << endl;
		}
	}
}

//------------------------------------------------------------------------------------------------------
//グラフの生成に関するルーチン
//------------------------------------------------------------------------------------------------------
//edge生成--------------------------------------------------------------------
void make_edge(GRAPH *g, OBSTACLES &obstacle){
	//エッジを保存するための領域確保
	int g_node_size = g->node.size();
	g->edge = vector<vector<EDGE> >(g_node_size);
	//エッジの追加
	for (int j = 0; j < g_node_size; j++){//From i To jにエッジがはれるか
		//ノードjの円内に中心を持つノードにエッジを張る．
		int xin = g->node[j].x;
		int yin = g->node[j].y;
		//エッジを列挙するために９近傍を参照
		int xnum_min, xnum_max;
		int ynum_min, ynum_max;
		g->calc_neighborhood(xin, yin, &xnum_min, &xnum_max, &ynum_min, &ynum_max);
		for (int k = ynum_min; k <= ynum_max; k++){
			for (int l = xnum_min; l <= xnum_max; l++){
				int indmax = g->data[k][l].size();
				for (int m = 0; m < indmax; m++){
					//九近傍内にあるノードの番号をiに代入
					int i = g->data[k][l][m];
					double dx = g->node[i].x - g->node[j].x;
					double dy = g->node[i].y - g->node[j].y;
					double nr = g->node[j].r * g->node[j].r;

					if (dx * dx + dy * dy < nr){//ノードjの円にノードiの中心があるなら
						int node_cost = pow(2, g->node[j].level);
						double v = g->node[j].v(g->node[i].x, g->node[i].y);
						double edge_cost = v + 1000.0;// / node_cost;
						g->edge[i].push_back(EDGE(j, edge_cost));
					}
				}
			}
		}
	}
}

//node生成--------------------------------------------------------------------
//壁と近いとき1を返す．また，obstaclex, obstacleyに近い壁の座標を出力する．
//出力される壁は一番最初に見つかったもので，最近近傍のものとは限らない．
int check_wall(const OBSTACLES *obstacle, double xin, double  yin, double rin){
	int xnum_min, xnum_max;
	int ynum_min, ynum_max;
	obstacle->calc_neighborhood(xin, yin, &xnum_min, &xnum_max, &ynum_min, &ynum_max);//9近傍を列挙

	for (int j = ynum_min; j <= ynum_max; j++){
		for (int i = xnum_min; i <= xnum_max; i++){
			for (int k = 0; k < obstacle->data[j][i].size(); k++){
				double obx = obstacle->data[j][i][k].x;
				double oby = obstacle->data[j][i][k].y;
				if (point_len(xin, obx, yin, oby) < rin * rin){
					return 1;
				}
			}
		}
	}
	return 0;
}
void make_node_sub(GRAPH *g, const OBSTACLES *obstacle, double xin, double yin, double r, int level){
	int iswall = check_wall(obstacle, xin, yin, r);
	if (iswall == 0){
		//障害物とかぶってなければノード追加
		g->add_node(NODE(xin, yin, r, level));
	}
	else{
		if (level < level_max){
			double r_next = r / sqrt(2.0) * 1.1;//1.1はマージン
			double d = r / 2.0;
			make_node_sub(g, obstacle, xin + d, yin + d, r_next, level + 1);
			make_node_sub(g, obstacle, xin - d, yin + d, r_next, level + 1);
			make_node_sub(g, obstacle, xin + d, yin - d, r_next, level + 1);
			make_node_sub(g, obstacle, xin - d, yin - d, r_next, level + 1);
		}
	}
}
void make_node(GRAPH *g, const OBSTACLES *obstacle){
	double d = discretization_width;//この幅で代表点を配置
	double r = r_first_level;//最初のレベルの円の半径
	int size_x = obstacle->xmax / d + 1.0;
	int size_y = obstacle->ymax / d + 1.0;

	for (int y = 1; y < size_y; y++){
		for (int x = 1; x < size_x; x++){
			double xin = discretization_width * x;
			double yin = discretization_width * y;
			//xin,yin,rの円を配置
			make_node_sub(g, obstacle, xin, yin, r, 0);
			if ((x - 1.0 / 2.0) < size_x || (y - 1.0 / 2.0) < size_y){
				xin = discretization_width * (x - 1.0 / 2.0);
				yin = discretization_width * (y - 1.0 / 2.0);
			}
			make_node_sub(g, obstacle, xin, yin, r, 0);
		}
	}
}
//------------------------------------------------------------------------------------------------------
//bellman_ford法で経路コスト計算
//------------------------------------------------------------------------------------------------------
bool bellman_ford(const int s, const GRAPH g, vector<double> &dist){//nは頂点数、sは開始頂点
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
//------------------------------------------------------------------------------------------------------
//Dijkstra法（優先度付きキュー実装）で経路コスト計算
//https://ja.wikipedia.org/wiki/%E3%83%80%E3%82%A4%E3%82%AF%E3%82%B9%E3%83%88%E3%83%A9%E6%B3%95#.E3.82.A2.E3.83.AB.E3.82.B4.E3.83.AA.E3.82.BA.E3.83.A0.E3.81.AE.E8.A7.A3.E8.AA.AC
//http://joints.blog111.fc2.com/blog-entry-28.html?sp を参考に．
//------------------------------------------------------------------------------------------------------
vector<int> dijkstra_withPriority(int s, const GRAPH g, vector<double> &dist)
{
	int g_node_size = g.node.size();
	dist = vector<double>(g_node_size, INF);//各ノードのコストを保存
	vector<int> prev(g_node_size);//手前のノードは何か？
	priority_queue<pair<int, int> > ensemble;//（ノード番号，コスト）をペアとした優先順位つきキュー

	dist[s] = 0;
	ensemble.push(pair<int, int>(s, 0));

	while (ensemble.size() != 0){
		int at = ensemble.top().first;
		ensemble.pop();
		for (int t = 0; t < g.edge[at].size(); t++){
			double edgecost = g.edge[at][t].cost;
			int i = g.edge[at][t].to;
			if (dist[i] > dist[at] + edgecost){
				dist[i] = dist[at] + edgecost;
				prev[i] = at;
				ensemble.push(pair<int, int>(i, dist[i]));
			}
		}
	}
	return prev;
}
//------------------------------------------------------------------------------------------------------
//制御則の運用に関するルーチン
//------------------------------------------------------------------------------------------------------
//レイヤー群の最小をとる----------------------------------------------------------
double v_min(GRAPH *g, double xin, double yin, vector<double> &dist, int *hoji_i){
	double v;//V(xin, yin)=min_{すべてのノードi}ノードiの関数V(xin,Yin)+始点sからノードiまでの重み
	double min = INF;
	//点xin,yinに重なっている可能性があるのは９近傍内の円のみ．
	//九近傍内の円を参照
	int xnum_min, xnum_max;
	int ynum_min, ynum_max;
	g->calc_neighborhood(xin, yin, &xnum_min, &xnum_max, &ynum_min, &ynum_max);
	for (int k = ynum_min; k <= ynum_max; k++){
		for (int l = xnum_min; l <= xnum_max; l++){
			int indmax = g->data[k][l].size();
			for (int m = 0; m < indmax; m++){
				//九近傍内にあるノードの番号をiに代入
				int p = g->data[k][l][m];
				//ノード番号pに対して演算
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
		}
	}
	return min;
}
//オイラー法による軌道計算オブジェクト---------------------------------------------
void dissasembled_differential(GRAPH *g, double x[2], double p[2], vector<double> &dist){
	int i;
	v_min(g, x[0], x[1], dist, &i);

	double xt = x[0] - g->node[i].x;
	double yt = x[1] - g->node[i].y;
	double sqr = sqrt(xt * xt + yt * yt);
	double k = 1.0 / sqr / pow(cos(M_PI / 2.0 / g->node[i].r), 2.0);

	p[0] = xt * k;
	p[1] = yt * k;
}

// pをuに変換---------------------------------------------------------------------
void euler(GRAPH *g, double x[2], double xdot[2], vector<double> &dist){
	double p[2];
	double u[2];
	//制御則の計算
	dissasembled_differential(g, x, p, dist);
	u[0] = -p[0];
	u[1] = -p[1];
	//ダイナミクスの計算
	xdot[0] = u[0];
	xdot[1] = u[1];
}

//オイラー法本体---------------------------------------------------------------------
void calc_euler(GRAPH *g, vector<double> &dist, int s){
	double log_next_time = 0;//間引きながらログを表示する用
	double x[2] = { start[0], start[1] };//状態量

	ofstream outputfile_xy(Fxy);

	for (double time = 0.0; time <= euler_end_time; time += DT){
		//ダイナミクスの更新
		double xdot[2];
		euler(g, x, xdot, dist);

		//間引きながらログ表示
		if (time > log_next_time){
			log_next_time += 1e-2;
			cout << x[0] << " " << x[1] << endl;
			outputfile_xy << x[0] << "," << x[1] << "," << time << endl;
		}

		//状態量の更新
		x[0] += DT * xdot[0];
		x[1] += DT * xdot[1];
	}

	outputfile_xy.close();
}

//Ｖのグラフデータを出力----------------------------------------------------------

void print_graph(GRAPH *g, const OBSTACLES &obstacle, vector<double> &dist){
	double xmin = obstacle.xmin;
	double xmax = obstacle.xmax;
	double ymin = obstacle.ymin;
	double ymax = obstacle.ymax;
	double vmax = 60000;
	double meshwidth_x = (xmax - xmin) / 100.0;
	double meshwidth_y = (ymax - ymin) / 100.0;

	FILE *fp;
	double x, y, z;
	fp = fopen("result_3d_z.txt", "w");
	for (y = ymin; y < ymax; y += meshwidth_y){
		for (x = xmin; x < xmax; x += meshwidth_x){
			int i;
			z = v_min(g, x, y, dist, &i);
			if (z > vmax)z = vmax;
			fprintf(fp, "%lf ", z);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	fp = fopen("result_3d_x.txt", "w");
	for (x = xmin; x < xmax; x += meshwidth_x){
		fprintf(fp, "%lf ", x);
	}
	fclose(fp);
	fp = fopen("result_3d_y.txt", "w");
	for (y = ymin; y < ymax; y += meshwidth_y){
		fprintf(fp, "%lf ", y);
	}
	fclose(fp);
}
//main--------------------------------------------------------------------------
int main() {
	read_simulation_condition();        //条件ファイルよりシミュレーション条件読み込み

	int s;//ゴールノード番号
	vector<double> dist; // 最短距離
	OBSTACLES obstacle(obstacle_file);
	GRAPH g(obstacle);

	cout << "makenode start" << endl;
	make_node(&g, &obstacle);
	cout << "makenode end" << endl;

	cout << "makeedge start" << endl;
	make_edge(&g, obstacle);                 //エッジ作成
	cout << "makeedge end" << endl;

	cout << "nodenum : " << g.node.size() << endl;
	int edgenum = 0;
	for (int i = 0; i < g.edge.size(); i++)edgenum += g.edge[i].size();
	cout << "edgenum : " << edgenum << endl;

	g.calc_nodenumber(goal[0], goal[1], &s); //与えられた座標に一番近い中心を持つノード番号を計算

	cout << "Graph algorith start" << endl;
	//	bellman_ford(s, g, dist);            //最短経路計算
	dijkstra_withPriority(s, g, dist);       //最短経路計算
	cout << "Graph algorith end" << endl;

	print_graph(&g, obstacle, dist);         //グラフデータ出力
	g.print();                               //離散点の２D画像ファイルに出力
	calc_euler(&g, dist, s);                 //移動経路計算
	return 0;
}
