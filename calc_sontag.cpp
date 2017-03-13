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
#include <queue>//Dijkstra�A���S���Y���̂��߂ɕK�v

using namespace std;
#define point_len(x, xx, y, yy)((x - xx) * (x - xx) + (y - yy) * (y - yy))//�_�Ɠ_�̋���
const double INF = DBL_MAX;             //double�^�̍ő�l

//------------------------------------------------------------------------------------------------------
//�V�~�����[�V���������Ɋւ���ݒ�
//------------------------------------------------------------------------------------------------------
const char *simulation_condition = "condition_lab.txt";//�V�~�����[�V�����������L�q���Ă���t�@�C����
//const char *simulation_condition = "condition_simple.txt";//�V�~�����[�V�����������L�q���Ă���t�@�C����

//------------------------------------------------------------------------------------------------------
//�I�C���[�@�Ɋւ���p�����[�^
//------------------------------------------------------------------------------------------------------
const char *Fxy = "resultXY.csv";       //���O�̏o�͐�
const double euler_end_time = 400.0;    //�V�~�����[�V�����̍ő厞��
const double DT = 0.0001;               //�I�C���[�@�̎��ԃp�����[�^
double start[2];                        //������ԗ�

//------------------------------------------------------------------------------------------------------
//����ڕW�Ɋւ���p�����[�^
//------------------------------------------------------------------------------------------------------
double goal[2];                         //�ڕW�ʒu

//------------------------------------------------------------------------------------------------------
//��Ԃ̗��U���Ɋւ���ݒ�
//------------------------------------------------------------------------------------------------------
char obstacle_file[256];                //��Q���f�[�^���L�^����Ă���e�L�X�g�t�@�C���̃t�@�C����
const char *Fname_matlab = "result.csv";//���U�_�Q��matlab�̃t�H�[�}�b�g�ŏo�͂���ۂ̃t�@�C����
const double discretization_width = 16;  //���U���̑�1���x���ɂ�����e��(�_���ǂ̂��炢�̕��Ŕz�u���邩)
const double r_first_level = 12.0;       //��1���x���̗��U���ɂ�����~�̔��a.
//r_first_level > sqrt(2.0)/2.0 * discretization_width�𖞂����K�v����
const int level_max = 3;                //���U���̍ċN�v���Z�X�����񂭂肩������

//------------------------------------------------------------------------------------------------------
//�����t�@�C�����V�~�����[�V����������ǂݍ���
//------------------------------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------------------------------
//��Q������舵�����߂̃N���X
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
	vector< vector< vector<xypoint> > > data;//�O�����z��ɏ�Q���̈ʒu����ۑ�
	//������
	OBSTACLES(const char *str);
	//�^����ꂽ�i���C���j���W����������̋�ߖT�i�O�㍶�E�΂߁j�̋����W��Ԃ��֐�
	//const�C���q��ǉ����Ă��āC�������茾���΁C�u�����ϐ���������ȁI�v�Ǝw�肵�Ă���D
	void calc_neighborhood(double xin, double yin, int *xoutmin, int *xoutmax, int *youtmin, int *youtmax) const;
	//����̃O���b�h�ɂ���_��S�ė񋓂���
	void print_obstacle(int y, int x) const;

};

OBSTACLES::OBSTACLES(const char *str){
	lwidth = r_first_level;
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
//const�C���q��ǉ����Ă��āC�������茾���΁C�u�����ϐ���������ȁI�v�Ǝw�肵�Ă���D
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
//����̃O���b�h�ɂ���_��S�ė񋓂���
void OBSTACLES::print_obstacle(int y, int x) const{
	for (int i = 0; i < data[y][x].size(); i++){
		cout << "obstacle: " << data[y][x][i].x << " " << data[y][x][i].y << endl;
	}
}
//------------------------------------------------------------------------------------------------------
//��Ė@�ɂ�����O���t�\������舵�����߂̃N���X
//------------------------------------------------------------------------------------------------------
//�m�[�h�̒�`-------------------------------------------------------------------
class NODE{
public:
	double x, y, r;//( x���W, y���W, ���a,
	int level;//�����x���ڂ̍ċN�A���S���Y���Œǉ����ꂽ��

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

//�G�b�W�̒�`-------------------------------------------------------------------
class EDGE{
public:
	int to;
	double cost;
	//�i�ǂ̃m�[�h�ɂȂ��邩, �G�b�W�R�X�g�j
	EDGE(int to_, double cost_) : to(to_), cost(cost_) {}
};

//�O���t�̒�`-------------------------------------------------------------------
class GRAPH{
public:
	//�O���t�f�[�^�\��
	vector<NODE> node;
	vector<vector<EDGE> > edge;
	//�N�C�b�N�A�N�Z�X�̂��߂̃f�[�^�\��
	double xmin, xmax, ymin, ymax;
	double lwidth;
	int xnum, ynum;
	vector< vector< vector<int> > > data;//�O�����z��Ƀm�[�h�̃C���f�b�N�X��ۑ�

	//������
	GRAPH(OBSTACLES &obstacle);
	//�m�[�h�ǉ�
	void add_node(NODE n);
	//�^����ꂽ�i���C���j���W����������̋�ߖT�i�O�㍶�E�΂߁j�̋����W��Ԃ��֐�
	//const�C���q��ǉ����Ă��āC�������茾���΁C�u�����ϐ���������ȁI�v�Ǝw�肵�Ă���D
	void calc_neighborhood(double xin, double yin, int *xoutmin, int *xoutmax, int *youtmin, int *youtmax) const;
	//�^����ꂽ���W�Ɉ�ԋ߂����S�����m�[�h�ԍ���Ԃ�
	void calc_nodenumber(const double x, const double y, int *s);
	void print(void);//�O���t�̏��\��
	void print_full(void);//�O���t�̏��\��
};
//������------------------------------------------------------------------------
GRAPH::GRAPH(OBSTACLES &obstacle){
	xmin = obstacle.xmin;
	xmax = obstacle.xmax;
	ymin = obstacle.ymin;
	ymax = obstacle.ymax;
	lwidth = obstacle.lwidth;
	//�s��̗v�f�����v�Z
	xnum = (int)((xmax - xmin) / lwidth) + 1;
	ynum = (int)((ymax - ymin) / lwidth) + 1;
	//�z��̃T�C�Y�ύX
	data.resize(ynum);		// ()���̐������v�f���ɂȂ�
	for (int i = 0; i < ynum; i++){
		data[i].resize(xnum);
	}
}
//�m�[�h�ǉ�--------------------------------------------------------------------
void GRAPH::add_node(NODE n){
	int current_index = node.size();//�ǉ�������̂����Ԗڂ�index��
	node.push_back(n);//�f�[�^�ɂo�t�r�g
	//�ǂ̋��ɓ��邩�v�Z
	int y = (int)((n.y - ymin) / lwidth);
	int x = (int)((n.x - xmin) / lwidth);
	//���Y���Ƀm�[�h��index��push
	if (0 <= y && y < ynum && 0 <= x && x < xnum){
		data[y][x].push_back(current_index);
	}
}
//�^����ꂽ�i���C���j���W����������̋�ߖT�i�O�㍶�E�΂߁j�̋����W��Ԃ��֐�
//const�C���q��ǉ����Ă��āC�������茾���΁C�u�����ϐ���������ȁI�v�Ǝw�肵�Ă���D
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
//�^����ꂽ���W�Ɉ�ԋ߂����S�����m�[�h�ԍ���Ԃ�--------------------------------
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
//���U�_�Q���o�́Dmatlab�Ő}���m�F�D----------------------------------------------
void GRAPH::print(void){
	ofstream outputfile(Fname_matlab);

	for (int i = 0; i < node.size(); i++){
		outputfile << node[i].x << ", " << node[i].y << ", " << node[i].r << ", " << node[i].level << endl;
	}
	outputfile.close();
}
//�O���t�̏��\��------------------------------------------------------------
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
//�O���t�̐����Ɋւ��郋�[�`��
//------------------------------------------------------------------------------------------------------
//edge����--------------------------------------------------------------------
void make_edge(GRAPH *g, OBSTACLES &obstacle){
	//�G�b�W��ۑ����邽�߂̗̈�m��
	int g_node_size = g->node.size();
	g->edge = vector<vector<EDGE> >(g_node_size);
	//�G�b�W�̒ǉ�
	for (int j = 0; j < g_node_size; j++){//From i To j�ɃG�b�W���͂�邩
		//�m�[�hj�̉~���ɒ��S�����m�[�h�ɃG�b�W�𒣂�D
		int xin = g->node[j].x;
		int yin = g->node[j].y;
		//�G�b�W��񋓂��邽�߂ɂX�ߖT���Q��
		int xnum_min, xnum_max;
		int ynum_min, ynum_max;
		g->calc_neighborhood(xin, yin, &xnum_min, &xnum_max, &ynum_min, &ynum_max);
		for (int k = ynum_min; k <= ynum_max; k++){
			for (int l = xnum_min; l <= xnum_max; l++){
				int indmax = g->data[k][l].size();
				for (int m = 0; m < indmax; m++){
					//��ߖT���ɂ���m�[�h�̔ԍ���i�ɑ��
					int i = g->data[k][l][m];
					double dx = g->node[i].x - g->node[j].x;
					double dy = g->node[i].y - g->node[j].y;
					double nr = g->node[j].r * g->node[j].r;

					if (dx * dx + dy * dy < nr){//�m�[�hj�̉~�Ƀm�[�hi�̒��S������Ȃ�
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

//node����--------------------------------------------------------------------
//�ǂƋ߂��Ƃ�1��Ԃ��D�܂��Cobstaclex, obstacley�ɋ߂��ǂ̍��W���o�͂���D
//�o�͂����ǂ͈�ԍŏ��Ɍ����������̂ŁC�ŋߋߖT�̂��̂Ƃ͌���Ȃ��D
int check_wall(const OBSTACLES *obstacle, double xin, double  yin, double rin){
	int xnum_min, xnum_max;
	int ynum_min, ynum_max;
	obstacle->calc_neighborhood(xin, yin, &xnum_min, &xnum_max, &ynum_min, &ynum_max);//9�ߖT���

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
		//��Q���Ƃ��Ԃ��ĂȂ���΃m�[�h�ǉ�
		g->add_node(NODE(xin, yin, r, level));
	}
	else{
		if (level < level_max){
			double r_next = r / sqrt(2.0) * 1.1;//1.1�̓}�[�W��
			double d = r / 2.0;
			make_node_sub(g, obstacle, xin + d, yin + d, r_next, level + 1);
			make_node_sub(g, obstacle, xin - d, yin + d, r_next, level + 1);
			make_node_sub(g, obstacle, xin + d, yin - d, r_next, level + 1);
			make_node_sub(g, obstacle, xin - d, yin - d, r_next, level + 1);
		}
	}
}
void make_node(GRAPH *g, const OBSTACLES *obstacle){
	double d = discretization_width;//���̕��ő�\�_��z�u
	double r = r_first_level;//�ŏ��̃��x���̉~�̔��a
	int size_x = obstacle->xmax / d + 1.0;
	int size_y = obstacle->ymax / d + 1.0;

	for (int y = 1; y < size_y; y++){
		for (int x = 1; x < size_x; x++){
			double xin = discretization_width * x;
			double yin = discretization_width * y;
			//xin,yin,r�̉~��z�u
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
//bellman_ford�@�Ōo�H�R�X�g�v�Z
//------------------------------------------------------------------------------------------------------
bool bellman_ford(const int s, const GRAPH g, vector<double> &dist){//n�͒��_���As�͊J�n���_
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
//------------------------------------------------------------------------------------------------------
//Dijkstra�@�i�D��x�t���L���[�����j�Ōo�H�R�X�g�v�Z
//https://ja.wikipedia.org/wiki/%E3%83%80%E3%82%A4%E3%82%AF%E3%82%B9%E3%83%88%E3%83%A9%E6%B3%95#.E3.82.A2.E3.83.AB.E3.82.B4.E3.83.AA.E3.82.BA.E3.83.A0.E3.81.AE.E8.A7.A3.E8.AA.AC
//http://joints.blog111.fc2.com/blog-entry-28.html?sp ���Q�l�ɁD
//------------------------------------------------------------------------------------------------------
vector<int> dijkstra_withPriority(int s, const GRAPH g, vector<double> &dist)
{
	int g_node_size = g.node.size();
	dist = vector<double>(g_node_size, INF);//�e�m�[�h�̃R�X�g��ۑ�
	vector<int> prev(g_node_size);//��O�̃m�[�h�͉����H
	priority_queue<pair<int, int> > ensemble;//�i�m�[�h�ԍ��C�R�X�g�j���y�A�Ƃ����D�揇�ʂ��L���[

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
//���䑥�̉^�p�Ɋւ��郋�[�`��
//------------------------------------------------------------------------------------------------------
//���C���[�Q�̍ŏ����Ƃ�----------------------------------------------------------
double v_min(GRAPH *g, double xin, double yin, vector<double> &dist, int *hoji_i){
	double v;//V(xin, yin)=min_{���ׂẴm�[�hi}�m�[�hi�̊֐�V(xin,Yin)+�n�_s����m�[�hi�܂ł̏d��
	double min = INF;
	//�_xin,yin�ɏd�Ȃ��Ă���\��������̂͂X�ߖT���̉~�̂݁D
	//��ߖT���̉~���Q��
	int xnum_min, xnum_max;
	int ynum_min, ynum_max;
	g->calc_neighborhood(xin, yin, &xnum_min, &xnum_max, &ynum_min, &ynum_max);
	for (int k = ynum_min; k <= ynum_max; k++){
		for (int l = xnum_min; l <= xnum_max; l++){
			int indmax = g->data[k][l].size();
			for (int m = 0; m < indmax; m++){
				//��ߖT���ɂ���m�[�h�̔ԍ���i�ɑ��
				int p = g->data[k][l][m];
				//�m�[�h�ԍ�p�ɑ΂��ĉ��Z
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
//�I�C���[�@�ɂ��O���v�Z�I�u�W�F�N�g---------------------------------------------
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

// p��u�ɕϊ�---------------------------------------------------------------------
void euler(GRAPH *g, double x[2], double xdot[2], vector<double> &dist){
	double p[2];
	double u[2];
	//���䑥�̌v�Z
	dissasembled_differential(g, x, p, dist);
	u[0] = -p[0];
	u[1] = -p[1];
	//�_�C�i�~�N�X�̌v�Z
	xdot[0] = u[0];
	xdot[1] = u[1];
}

//�I�C���[�@�{��---------------------------------------------------------------------
void calc_euler(GRAPH *g, vector<double> &dist, int s){
	double log_next_time = 0;//�Ԉ����Ȃ��烍�O��\������p
	double x[2] = { start[0], start[1] };//��ԗ�

	ofstream outputfile_xy(Fxy);

	for (double time = 0.0; time <= euler_end_time; time += DT){
		//�_�C�i�~�N�X�̍X�V
		double xdot[2];
		euler(g, x, xdot, dist);

		//�Ԉ����Ȃ��烍�O�\��
		if (time > log_next_time){
			log_next_time += 1e-2;
			cout << x[0] << " " << x[1] << endl;
			outputfile_xy << x[0] << "," << x[1] << "," << time << endl;
		}

		//��ԗʂ̍X�V
		x[0] += DT * xdot[0];
		x[1] += DT * xdot[1];
	}

	outputfile_xy.close();
}

//�u�̃O���t�f�[�^���o��----------------------------------------------------------

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
	read_simulation_condition();        //�����t�@�C�����V�~�����[�V���������ǂݍ���

	int s;//�S�[���m�[�h�ԍ�
	vector<double> dist; // �ŒZ����
	OBSTACLES obstacle(obstacle_file);
	GRAPH g(obstacle);

	cout << "makenode start" << endl;
	make_node(&g, &obstacle);
	cout << "makenode end" << endl;

	cout << "makeedge start" << endl;
	make_edge(&g, obstacle);                 //�G�b�W�쐬
	cout << "makeedge end" << endl;

	cout << "nodenum : " << g.node.size() << endl;
	int edgenum = 0;
	for (int i = 0; i < g.edge.size(); i++)edgenum += g.edge[i].size();
	cout << "edgenum : " << edgenum << endl;

	g.calc_nodenumber(goal[0], goal[1], &s); //�^����ꂽ���W�Ɉ�ԋ߂����S�����m�[�h�ԍ����v�Z

	cout << "Graph algorith start" << endl;
	//	bellman_ford(s, g, dist);            //�ŒZ�o�H�v�Z
	dijkstra_withPriority(s, g, dist);       //�ŒZ�o�H�v�Z
	cout << "Graph algorith end" << endl;

	print_graph(&g, obstacle, dist);         //�O���t�f�[�^�o��
	g.print();                               //���U�_�̂QD�摜�t�@�C���ɏo��
	calc_euler(&g, dist, s);                 //�ړ��o�H�v�Z
	return 0;
}
