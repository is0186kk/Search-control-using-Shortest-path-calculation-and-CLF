//離散点の作成　障害物座標 複数
//四角く区切って効率化(壁との距離比較関数作成)
//計算時間算出
//
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#define map_size_x 50.0 //マップサイズ（横）
#define map_size_y 50.0 //マップサイズ（縦）
#define detail 4        //離散化の粗さを決める
#define density 4       //離散化の細かさを決める
#define rin 3.0         //全体離散点のノードの半径
#define F "test_dots.txt"
//座標(x, y)と座標(xx, yy)の距離を計算
#define point_len(x, xx, y, yy)((x - xx) * (x - xx) + (y - yy) * (y - yy))
using namespace std;

class DEST{//離散点群を格納
public:
	double dx, dy, dr;
	DEST(double d_x_, double d_y_, double d_r_) : dx(d_x_), dy(d_y_), dr(d_r_){}
};

class STAK{//ノード内に障害物を含む場合，この配列に格納．細かな離散化を行う．
public:
	double sx, sy, sr;
  STAK(double s_x_, double s_y_, double s_r_) : sx(s_x_), sy(s_y_), sr(s_r_){}
};

class MAP{//MAPクラス
public:
  vector<DEST> dest;
  vector<STAK> stak;

//離散点群を出力．matlabで図を確認．
	void print(){
		FILE *fp;
		fp = fopen("result.csv","w");
		for(int i = 0; i < dest.size(); i++){
			fprintf(fp, "%lf, %lf\n", dest[i].dx, dest[i].dy);//, dest[i].dr);
		}
		fclose(fp);
	}
};

class OBSTACLES{
private:
	typedef struct{
		double x;
		double y;
	}xypoint;
public:
	double xmin, xmax, ymin, ymax, rmax;
	vector< vector< vector<xypoint> > > data;//三次元配列に障害物の位置情報を保存
	//特定のグリッドにある点を全て列挙する
	void print_obstacle(int y, int x){
		for(int i = 0; i < data[y][x].size(); i++){
			cout << "obstacle: "<< data[y][x][i].x << " "<<data[y][x][i].y<< endl;
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
class timercounter{
	clock_t s;
	clock_t e;
public:
	void start(){
		s = clock();
	}
	void end(){
		e = clock();
	}
	void print(){
		double t = (double)(e - s) / CLOCKS_PER_SEC;
	        std::cout << "duration = " << t << "sec." << std::endl;
	}
	void end_print(){
		end();
		print();
	}
};

//eraseを使うための処理
template<typename T>
void remove(std::vector<T>& vector, unsigned int index){
	vector.erase(vector.begin() + index);
}
//壁に近すぎる離散点を消す
void erase_wall(MAP *map, OBSTACLES obstacle){
	for(int i = 0; i < map->dest.size(); i++){
		int y = (int)((map->dest[i].dy - obstacle.ymin) / obstacle.rmax);
		int x = (int)((map->dest[i].dx - obstacle.xmin) / obstacle.rmax);

		for(int w = 0; w < obstacle.data[y][x].size(); w++){
			if(point_len(map->dest[i].dx, obstacle.data[y][x][w].x,
				           map->dest[i].dy, obstacle.data[y][x][w].y ) < 0.001
								|| map->dest[i].dx < 0.0
								|| map->dest[i].dy < 0.0){
									 remove(map->dest, i);
									 i--;
			}
		}
	}
}

// 配列内に同じ座標が存在する場合にpush_backしないための関数
// 同じ座標に離散点を作ろうとしていた場合は半径が大きいほうを削除
int check_same(MAP *map, double x, double y, double r){
	for(int i = 0; i < map->dest.size(); i++){
		if(map->dest[i].dx == x && map->dest[i].dy == y){
			if(map->dest[i].dr > r){
				 map->dest[i].dr = r;
				return 0;
			}
			else return 1;
		}
	}
	return 0;
}

//障害物近傍の大きいノードを離散化
void make_small_node(MAP *map, double s_x, double s_y, double s_r, double wall_x, double wall_y){
  double s_xx, s_yy, s_rr;
  double r = log2(map->dest[0].dr / s_r);
  const double dx[5] = {0, -1, +1, -1, +1};
  const double dy[5] = {0, +1, +1, -1, -1};
  const double dz[4][2] = {{+1, 0}, {-1, 0}, {0, +1}, {0, -1}};
  double dz_x, dz_y;

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
        && point_len(dz_x, wall_x, dz_y, wall_y) <  s_r * s_r
        && point_len(dz_x, wall_x, dz_y, wall_y) > s_rr * s_rr
        && r < density - 1){
					if(check_same(map, dz_x, dz_y, s_rr) == 0)
					map->dest.push_back(DEST(dz_x, dz_y, s_rr));
				}
    }
		if(s_rr < map->dest[0].dr * pow(0.5, density)) break;//十分な細かさが実現したら終了
    else if(point_len(s_xx, wall_x, s_yy, wall_y) < s_rr * s_rr)//障害物をノード内に含むとき,より細かく離散化．
      make_small_node(map, s_xx, s_yy, s_rr, wall_x, wall_y);
    else{//障害物をノード内に含まないため,離散点とする．
			if(check_same(map, s_xx, s_yy, s_rr) == 0)
			map->dest.push_back(DEST(s_xx, s_yy, s_rr));
		}
  }
}

//壁と近いとき1を返す
int check_wall(MAP* map, OBSTACLES obstacle, double xin, double yin, double rrin, int a){
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
					else make_small_node(map, xin, yin, rrin,
															 obstacle.data[j][i][k].x, obstacle.data[j][i][k].y);
				}
			}
		}
	}
	return 0;
}

//平面全体の離散化
void desti_all2(MAP *map, OBSTACLES obstacle, double xin, double yin){
	if(check_wall(map, obstacle, xin, yin, rin, 0) == 0){
		if(check_same(map, xin, yin, rin) == 0)
		map->dest.push_back(DEST(xin, yin, rin));
		}
		else map->stak.push_back(STAK(xin, yin, rin));
}

//平面全体の離散化
void desti_all(MAP *map, OBSTACLES obstacle){
	double size_x = map_size_x / detail + 1.0;
	double size_y = map_size_y / detail + 1.0;
	for(int y = 1; y < size_y; y++){
		for(int x = 1; x < size_x; x++){
			double xin = detail * x;
			double yin = detail * y;

			desti_all2(map, obstacle, xin, yin);
			if((x - 1.0 / 2.0) < size_x  || (y - 1.0 / 2.0) < size_y){
			xin = detail * (x - 1.0 / 2.0);
			yin = detail * (y - 1.0 / 2.0);
			desti_all2(map, obstacle, xin, yin);
		}
		}
	}
}

//離散化を管理する関数
void discretization(MAP *map, OBSTACLES obstacle){
	desti_all(map, obstacle);
	for(int s = 0; s < map->stak.size(); s++){
		check_wall(map, obstacle, map->stak[s].sx, map->stak[s].sy, map->stak[s].sr, 1);
	}
	erase_wall(map, obstacle);
}
//データインポート
// void import(MAP *map){
//
//
//   for(int i = 0; i < n; i++){
// 		cin >>wall_x >>wall_y;
// 		map->wall.push_back(WALL(wall_x, wall_y));
// 	}
// }

int main(){
  MAP map;
	timercounter timer;

	double xmin = 0;
	double xmax = map_size_x;
	double ymin = 0;
	double ymax = map_size_y;
	double rmax = 10;

  //import(&map);
	double x, y, wall_x, wall_y, gx, gy;
	//int n;
	//初期値
  //cin >>x >>y;
	//cin >>n;
	//壁座標インポート
	OBSTACLES obstacle("test_dots.txt", xmin, xmax, ymin, ymax, rmax);
	//obstacle.print_obstacle(5, 5);//対象のグリッドに入る障害物一覧表示
	//cin >>gx >>gy;
	timer.start();
	discretization(&map, obstacle);
	timer.end_print();
	map.print();
}
