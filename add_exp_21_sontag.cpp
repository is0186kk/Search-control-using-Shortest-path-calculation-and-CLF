/*
���U���{sontag

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
#define density 4       //���U���ׂ̍��������߂�
#define rin 3.0         //�S�̗��U�_�̃m�[�h�̔��a
#define DT 0.01

using namespace std;

const int INF = 60000;

//class-------------------------------------------------------------------------
class NODE{
	public:
	double x, y, r;
	int a, b, c;
	NODE(double x_, double y_, double r_, int a_, int b_, int c_) : x(x_), y(y_), r(r_), a(a_), b(b_), c(c_){}//( x���W, y���W, ���a, ���U�_�ԍ�, �Ǎ��W���f��̍X�V����ɗp����)
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
	//�i�ǂ̃m�[�h�ɂȂ��邩, �G�b�W�R�X�g�j
	EDGE(int to_, double cost_) : to(to_), cost(cost_) {}
};

//�m�[�h���ɏ�Q�����܂ޏꍇ�C���̔z��Ɋi�[�D�ׂ��ȗ��U�����s����_�Ƃ���D
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
	void print(void);//�O���t�̏��\��
};

//���U�_�Q���o�́Dmatlab�Ő}���m�F�D----------------------------------------------
void GRAPH::print(void){
	ofstream outputfile(F);

	for(int i=0; i < node.size(); i++){
		outputfile<< node[i].x<<", "<< node[i].y<<", "<< node[i].r<<endl;
	}
	outputfile.close();
}
// //�O���t�̏��\��------------------------------------------------------------
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

//bellman_ford�@�Ōo�H�R�X�g�v�Z--------------------------------------------------
bool bellman_ford(int s, GRAPH g, vector<double> &dist){//n�͒��_���As�͊J�n���_
  dist = vector<double>(g.node.size(), INF);
  dist[s] = 0; // �J�n�_�̋�����0
  for (int i = 0; i < g.node.size(); i++) {
    for (int v = 0; v < g.node.size(); v++) {
      for (int k = 0; k < g.edge[v].size(); k++) {
			  EDGE e = g.edge[v][k];
        if (dist[v] != INF && dist[e.to] > dist[v] + e.cost){
          dist[e.to] = dist[v] + e.cost;
          if (i == g.node.size() - 1) return true; //n��ڂɂ��X�V������Ȃ畉�̕H������
        }
      }
	  }
  }
  return false;
}

//�_�Ɠ_�̋���-------------------------------------------------------------------
double point_len(double x, double xx, double y, double yy){
	double min_p_p = 500000.0;
	double length;

	return length = pow( (x - xx) * (x - xx) + (y - yy) * (y - yy), 0.5 );
}

//���C���[�Q�̍ŏ����Ƃ�
double v_min(GRAPH g, double xin, double yin, vector<double> &dist, int hoji_i[1]){
  double v = 0.0;//V(xin, yin)=min_{���ׂẴm�[�hi}�m�[�hi�̊֐�V(xin,Yin)+�n�_s����m�[�hi�܂ł̏d��
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

//�I�C���[�@�ɂ��O���v�Z�I�u�W�F�N�g---------------------------------------------
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

// p��u�ɕϊ�---------------------------------------------------------------------
void euler(GRAPH g, double x[2], double u[2], vector<double> &dist){
	double p[2];
	dissasembled_differential(g, x, p, dist);
	u[0] = -p[0];
	u[1] = -p[1];
}

//�I�C���[�@---------------------------------------------------------------------
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

//graph�`��t�@�C������----------------------------------------------------------
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
					&& g.wall[k].wy - 0.001 <= y && y <= g.wall[k].wy + 0.001) z = INF;//�Ȃ��Ă�������
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

//edge����----------------------------------------------------------------------
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
//�֐�remove�̒�`---------------------------------------------------------------
template<typename T>//remove(g->node, i)�@g->node��i�Ԗڂ��폜
void remove(std::vector<T>& vector, unsigned int index){
	vector.erase(vector.begin() + index);
}
//�w�肵�����U�_���폜-----------------------------------------------------------
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
// //�������ǂɋ߂����闣�U�_������΁A�����ō폜�@�@�Ȃ��Ă�������
// void erase_close_node(GRAPH *g){
// 	for(int i = 0; i < g->node.size(); i++){
// 		for(int j = 0; j < g->wall.size(); j++){
// 			if(point_len(g->node[i].x, g->wall[j].wx, g->node[i].y, g->wall[j].wy) < g->node[i].r){
// 				remove(g->node, i);
// 			}
// 		}
// 	}
// }
// �z����ɓ������W�����݂���ꍇ��push_back���Ȃ����߂̊֐�------------------------
// �������W�ɗ��U�_����낤�Ƃ��Ă����ꍇ�͔��a���傫���ق����폜
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

//��Q���ߖT�̑傫���m�[�h�𗣎U��-------------------------------------------------
void make_small_node(GRAPH *g, double s_x, double s_y, double s_r, double wall_x, double wall_y){
	double s_xx, s_yy, s_rr;
	int ain = 0, bin = -1, ccin = 0;
  double r = log2(g->node[0].r / s_r);
  const double dx[5] = {0, -1, +1, -1, +1};
  const double dy[5] = {0, +1, +1, -1, -1};
  const double dz[4][2] = {{+1, 0}, {-1, 0}, {0, +1}, {0, -1}};
  double dz_x, dz_y;
	ain = r;
	//���S�C����C�E��C�����C�E���@�̏��ɑ傫���~�̏�ɏ������~����ׂ�D
	for(int i = 0; i < 5; i++){
		s_rr = s_r * 0.5;
    s_xx = s_x + dx[i] * pow(0.5, r);
    s_yy = s_y + dy[i] * pow(0.5, r);

  //��Q�������ɍ���_�Q,�E��_�Q,�����_�Q,�E���_�Q���ł���̂�,�Q���m���q�����U�_��݌v�@
  //���U�����ׂ�������,��Q���ɂ��߂����Ȃ��Ƃ�,�������W�����苗���̓_�𗣎U�_�Ƃ���D
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
			if(s_rr < g->node[0].r * pow(0.5, density)) break;//�\���ȍׂ���������������I��
			else if(point_len(s_xx, wall_x, s_yy, wall_y) < s_rr)//��Q�����m�[�h���Ɋ܂ނƂ�,���ׂ������U���D
			make_small_node(g, s_xx, s_yy, s_rr, wall_x, wall_y);
			else{//��Q�����m�[�h���Ɋ܂܂Ȃ����ߗ��U�_�Ƃ���D
				if(check_same(g, s_xx, s_yy, s_rr) == 0)
				g->node.push_back(NODE(s_xx, s_yy, s_rr, ain, bin, ccin));
			}
		}
	}


//���ʑS�̗̂��U��---------------------------------------------------------------
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

//���U�����Ǘ�����֐�-----------------------------------------------------------
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

//�f�[�^�C���|�[�g--�����n����----------------------------------------------------
void import(GRAPH *g, double start[2]){
  double xin, yin, w_x, w_y;
	int ain = 0, bin = -1, ccin = 0;
	int N;

	cin >> xin >> yin;
	cin >> N;
	start[0] = xin;
	start[1] = yin;

	//�ǈ���Ȃ��Ƃ��́C(-100, -100)�i���{�b�g���s�����Ȃ����W�j�ɕǂ�����Ƃ݂Ȃ��D
	//g.wall.size()�̃��[�v�����邽��
	if(N == 0)g->wall.push_back(WALL(-100.0, -100.0));

	//�Ǎ��W�C���|�[�g
	for(int n = 0; n < N; n++){
		cin >> w_x >> w_y;
		g->wall.push_back(WALL(w_x, w_y));
	}
}
//�_�Ɠ_�̋����ŃS�[���m�[�h�ԍ����擾---------------------------------------------
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
  vector<double> dist; // �ŒZ����

  GRAPH g;
  import(&g, start);//�t�@�C���������
	discretization(&g);//��ԗ��U��
	goal_import(&g, &s);//�ڕW���W�ݒ�
	g.edge=vector<vector<EDGE> >(g.node.size());
	make_edge(&g);//�G�b�W�쐬
	bellman_ford(s, g, dist);//�ŒZ�o�H���o
	//graph(g, dist);//�RD�摜�t�@�C���ɏo��
  g.print();//���U�_�̂QD�摜�t�@�C���ɏo��
	calc_euler(g, start, dist, s);//�ړ��o�H�v�Z
  return 0;
}
