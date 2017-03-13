/*
���U���{sontag

���̓t�H�[�}�b�g

�������V�~�����[�V�����ݒ���L�q����e�L�X�g�t�@�C���̓��e������
��Q�������L�q����e�L�X�g�t�@�C���̃t�@�C����
���{�b�g�̏����ʒu�̂����W�@�����W
�ڕW�ʒu�̂����W�@�����W

�������V�~�����[�V�����ݒ���L�q����e�L�X�g�t�@�C���̓��e������

��������Q�������L�q����e�L�X�g�t�@�C���̓��e������
map_size_xmin map_size_ymin
map_size_xmax map_size_ymax
��Q��1�̂����W�@��Q��1�̂����W
��Q��2�̂����W�@��Q��2�̂����W
�E
�E
�E
��������Q�������L�q����e�L�X�g�t�@�C���̓��e�����܂Ł�����


*/
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <float.h>//double�^�̍ő�l���ۑ�����Ă���

using namespace std;
#define point_len(x, xx, y, yy)((x - xx) * (x - xx) + (y - yy) * (y - yy))//�_�Ɠ_�̋���
const double INF = DBL_MAX;             //double�^�̍ő�l

//------------------------------------------------------------------------------------------------------
//�V�~�����[�V���������Ɋւ���ݒ�
//------------------------------------------------------------------------------------------------------
//const char *simulation_condition = "condition_lab.txt";//�V�~�����[�V�����������L�q���Ă���t�@�C����
const char *simulation_condition = "condition_simple.txt";//�V�~�����[�V�����������L�q���Ă���t�@�C����

//------------------------------------------------------------------------------------------------------
//��Ԃ̗��U���Ɋւ���ݒ�
//------------------------------------------------------------------------------------------------------
char obstacle_file[256];               //��Q���f�[�^���L�^����Ă���e�L�X�g�t�@�C���̃t�@�C����
const char *Fname_matlab = "result.csv";//���U�_�Q��matlab�̃t�H�[�}�b�g�ŏo�͂���ۂ̃t�@�C����
const double detail = 4;                //���U���̑�1���x���ɂ�����e��(�_���ǂ̂��炢�̕��Ŕz�u���邩)
const double rin = 3.0;                 //��1���x���̗��U���ɂ�����~�̔��a.rin > sqrt(2.0)/2.0 * detail�𖞂����K�v����
const int density = 2;                  //���U���̍ċN�v���Z�X�����񂭂肩������
const double K = (int)(rin + 0.5);        //K>rin�ł���K�v������
//------------------------------------------------------------------------------------------------------
//�I�C���[�@�Ɋւ���p�����[�^
//------------------------------------------------------------------------------------------------------
const char *Fxy = "resultXY.csv";       //���O�̏o�͐�
const double euler_end_time = 100.0;    //�V�~�����[�V�����̍ő厞��
const double DT = 0.01;                 //�I�C���[�@�̎��ԃp�����[�^
double start[2];                        //������ԗ�
//------------------------------------------------------------------------------------------------------
//����ڕW�Ɋւ���p�����[�^
//------------------------------------------------------------------------------------------------------
double goal[2];                         //�ڕW�ʒu

//class-------------------------------------------------------------------------
//a : ���U�����ׂ������邩�̕~���l�C�O�͑S�̗̂��U���i���x���P�j�@�P�̓��x���Q�ȏ�
//b :
//c :
class NODE{
public:
	double x, y, r;
	int a, b, c;
	NODE(double x_, double y_, double r_, int a_, int b_) : x(x_), y(y_), r(r_), a(a_), b(b_){}//( x���W, y���W, ���a, ���U�_�ԍ�, �Ǎ��W���f��̍X�V����ɗp����)
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
	//�i�ǂ̃m�[�h�ɂȂ��邩, �G�b�W�R�X�g�j
	EDGE(int to_, double cost_) : to(to_), cost(cost_) {}
};


class GRAPH{
public:
	vector<NODE> node;
	vector<vector<EDGE> > edge;
	//�^����ꂽ���W�Ɉ�ԋ߂����S�����m�[�h�ԍ���Ԃ�
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

	void print(void);//�O���t�̏��\��
};

//���U�_�Q���o�́Dmatlab�Ő}���m�F�D----------------------------------------------
void GRAPH::print(void){
	ofstream outputfile(Fname_matlab);

	for (int i = 0; i < node.size(); i++){
		outputfile << node[i].x << ", " << node[i].y << ", " << node[i].r << ", " << node[i].a << endl;
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


class OBSTACLES{
private:
	typedef struct{
		double x;
		double y;
	}xypoint;
public:
	double xmin, xmax, ymin, ymax;
	double lwidth;
	vector< vector< vector<xypoint> > > data;//�O�����z��ɏ�Q���̈ʒu����ۑ�

	//����̃O���b�h�ɂ���_��S�ė񋓂���
	void print_obstacle(int y, int x){
		for (int i = 0; i < data[y][x].size(); i++){
			cout << "obstacle: " << data[y][x][i].x << " " << data[y][x][i].y << endl;
		}
	}

	//������
	OBSTACLES(const char *str){
		lwidth = K;
		//�f�[�^�ǂݍ��ݏ���
		std::ifstream ifs(str);
		if (!ifs) {
			cerr << "file failure" << endl;
		}
		//�w�b�_�ǂݍ���
		ifs >> xmin >> xmax;
		ifs >> ymin >> ymax;

		//�f�[�^�{�̓ǂݍ��ݏ���
		std::cout << "Read obstacle data from " << str << std::endl;
		std::cout << "x coordinate " << xmin << "--" << xmax << std::endl;
		std::cout << "y coordinate " << ymin << "--" << ymax << std::endl;
		std::cout << "koushi haba " << lwidth << std::endl;
		//�s��̗v�f�����v�Z
		int xnum = (int)((xmax - xmin) / lwidth) + 1;
		int ynum = (int)((ymax - ymin) / lwidth) + 1;
		//�z��̃T�C�Y�ύX
		data.resize(ynum);		// ()���̐������v�f���ɂȂ�
		for (int i = 0; i < ynum; i++){
			data[i].resize(xnum);
		}

		//�f�[�^�{�̓ǂݍ���
		while (!ifs.eof()){
			double xin, yin;
			//�Ƃ肠��������
			ifs >> xin >> yin;
			//�ǂ̋��ɓ��邩�v�Z
			int y = (int)((yin - ymin) / lwidth);
			int x = (int)((xin - xmin) / lwidth);
			//���Y����push
			if (0 <= y && y < ynum && 0 <= x && x < xnum){
				xypoint inpoint = { xin, yin };
				data[y][x].push_back(inpoint);
			}
		}
	}
	//�^����ꂽ�i���C���j���W����������̋�ߖT�i�O�㍶�E�΂߁j�̋����W��Ԃ��֐�
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


//bellman_ford�@�Ōo�H�R�X�g�v�Z--------------------------------------------------
bool bellman_ford(int s, GRAPH g, vector<double> &dist){//n�͒��_���As�͊J�n���_
	int g_node_size = g.node.size();
	dist = vector<double>(g_node_size, INF);
	dist[s] = 0; // �J�n�_�̋�����0
	for (int i = 0; i < g_node_size; i++) {
		for (int v = 0; v < g_node_size; v++) {
			int g_edgev_size = g.edge[v].size();
			for (int k = 0; k < g_edgev_size; k++) {
				EDGE e = g.edge[v][k];
				if (dist[v] != INF && dist[e.to] > dist[v] + e.cost){
					dist[e.to] = dist[v] + e.cost;
					if (i == g_node_size - 1) return true; //n��ڂɂ��X�V������Ȃ畉�̕H������
				}
			}
		}
	}
	return false;
}

//�ǂƋ߂��Ƃ�1��Ԃ��D�܂��Cobstaclex, obstacley�ɋ߂��ǂ̍��W���o�͂���D
//�o�͂����ǂ͈�ԍŏ��Ɍ����������̂ŁC�ŋߋߖT�̂��̂Ƃ͌���Ȃ��D
int check_wall(GRAPH* g, OBSTACLES obstacle, double xin, double  yin, double rrin){
	int xnum_min, xnum_max;
	int ynum_min, ynum_max;
	obstacle.calc_neighborhood(xin, yin, &xnum_min, &xnum_max, &ynum_min, &ynum_max);//9�ߖT���

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

//���C���[�Q�̍ŏ����Ƃ�----------------------------------------------------------
double v_min(GRAPH *g, OBSTACLES obstacle,
	double xin, double yin, vector<double> &dist, int *hoji_i){
	double v;//V(xin, yin)=min_{���ׂẴm�[�hi}�m�[�hi�̊֐�V(xin,Yin)+�n�_s����m�[�hi�܂ł̏d��
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
//�I�C���[�@�ɂ��O���v�Z�I�u�W�F�N�g---------------------------------------------
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

// p��u�ɕϊ�---------------------------------------------------------------------
void euler(GRAPH *g, OBSTACLES obstacle, double x[2], double u[2], vector<double> &dist){
	double p[2];
	dissasembled_differential(g, obstacle, x, p, dist);
	u[0] = -p[0];
	u[1] = -p[1];
}

//�I�C���[�@�{��---------------------------------------------------------------------
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




//edge����--------------------------------------------------------------------
void make_edge(GRAPH *g, OBSTACLES obstacle){
	int g_node_size = g->node.size();
	g->edge = vector<vector<EDGE> >(g_node_size);
	for (int i = 0; i < g_node_size; i++){//From i To j�ɃG�b�W���͂�邩
		for (int j = 0; j < g_node_size; j++){
			double dx = g->node[i].x - g->node[j].x;
			double dy = g->node[i].y - g->node[j].y;
			double nr = g->node[j].r * g->node[j].r;

			if (dx * dx + dy * dy < nr){//�m�[�hj�̉~�Ƀm�[�hi�̒��S������Ȃ�
				int node_cost = pow(2, g->node[j].a);
				double v = g->node[j].v(g->node[i].x, g->node[i].y);
				double edge_cost = (v + 1000.0 / node_cost);
				g->edge[i].push_back(EDGE(j, edge_cost));
			}
		}
	}
}

//�����t�@�C�����V�~�����[�V����������ǂݍ���------------------------------------
void read_simulation_condition(void){
	std::ifstream ifs(simulation_condition);
	cout << simulation_condition << endl;
	if (!ifs) {
		cerr << "simulation file failure" << endl;
	}
	ifs.getline(obstacle_file, 256 - 1);//��Q���f�[�^
	ifs >> start[0] >> start[1];//�����ʒu
	ifs >> goal[0] >> goal[1];//�ڕW�ʒu
	cout << obstacle_file << endl;
}
void make_node_sub(GRAPH *g, OBSTACLES obstacle, double xin, double yin, double r, int level){
	int iswall = check_wall(g, obstacle, xin, yin, r);
	if (iswall == 0){
		//��Q���Ƃ��Ԃ��ĂȂ���΃m�[�h�ǉ�
		g->node.push_back(NODE(xin, yin, r, level, 0));
	}
	else{
		if (level < density){
			double r_next = r / sqrt(2.0) * 1.1;//1.1�̓}�[�W��
			double d = r / 2.0;
			make_node_sub(g, obstacle, xin + d, yin + d, r_next, level + 1);
			make_node_sub(g, obstacle, xin - d, yin + d, r_next, level + 1);
			make_node_sub(g, obstacle, xin + d, yin - d, r_next, level + 1);
			make_node_sub(g, obstacle, xin - d, yin - d, r_next, level + 1);
		}
	}
}
void make_node(GRAPH *g, OBSTACLES obstacle){
	double d = detail;//���̕��ő�\�_��z�u
	double r = d * 3.0 / 4.0;//�ŏ��̃��x���̉~�̔��a
	int size_x = obstacle.xmax / d + 1.0;
	int size_y = obstacle.ymax / d + 1.0;

	for (int y = 1; y < size_y; y++){
		for (int x = 1; x < size_x; x++){
			double xin = detail * x;
			double yin = detail * y;
			//xin,yin,r�̉~��z�u
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
	read_simulation_condition();        //�����t�@�C�����V�~�����[�V���������ǂݍ���

	int s;//�S�[���m�[�h�ԍ�
	vector<double> dist; // �ŒZ����
	OBSTACLES obstacle(obstacle_file);
	GRAPH g;

	cout << "makenode start" << endl;
	make_node(&g, obstacle);
	//	discretization(&g, obstacle);       //��ԗ��U���ɂ��m�[�h�쐬
	cout << "nodenum : " << g.node.size() << endl;
	cout << "makenode end" << endl;

	cout << "makeedge start" << endl;
	make_edge(&g, obstacle);            //�G�b�W�쐬
	cout << "makeedge end" << endl;

	g.calc_nodenumber(goal[0], goal[1], &s);  //�^����ꂽ���W�Ɉ�ԋ߂����S�����m�[�h�ԍ����v�Z

	cout << "bellmanford start" << endl;
	bellman_ford(s, g, dist);                 //�ŒZ�o�H�v�Z
	cout << "bellmanford end" << endl;

	g.print();                                //���U�_�̂QD�摜�t�@�C���ɏo��
	calc_euler(&g, obstacle, dist, s);        //�ړ��o�H�v�Z
	return 0;
}
