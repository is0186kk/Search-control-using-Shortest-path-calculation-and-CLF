//���U�_�̍쐬�@��Q�����W ����
//�l�p����؂��Č�����(�ǂƂ̋�����r�֐��쐬)
//�v�Z���ԎZ�o
//
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#define map_size_x 50.0 //�}�b�v�T�C�Y�i���j
#define map_size_y 50.0 //�}�b�v�T�C�Y�i�c�j
#define detail 4        //���U���̑e�������߂�
#define density 4       //���U���ׂ̍��������߂�
#define rin 3.0         //�S�̗��U�_�̃m�[�h�̔��a
#define F "test_dots.txt"
//���W(x, y)�ƍ��W(xx, yy)�̋������v�Z
#define point_len(x, xx, y, yy)((x - xx) * (x - xx) + (y - yy) * (y - yy))
using namespace std;

class DEST{//���U�_�Q���i�[
public:
	double dx, dy, dr;
	DEST(double d_x_, double d_y_, double d_r_) : dx(d_x_), dy(d_y_), dr(d_r_){}
};

class STAK{//�m�[�h���ɏ�Q�����܂ޏꍇ�C���̔z��Ɋi�[�D�ׂ��ȗ��U�����s���D
public:
	double sx, sy, sr;
  STAK(double s_x_, double s_y_, double s_r_) : sx(s_x_), sy(s_y_), sr(s_r_){}
};

class MAP{//MAP�N���X
public:
  vector<DEST> dest;
  vector<STAK> stak;

//���U�_�Q���o�́Dmatlab�Ő}���m�F�D
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
	vector< vector< vector<xypoint> > > data;//�O�����z��ɏ�Q���̈ʒu����ۑ�
	//����̃O���b�h�ɂ���_��S�ė񋓂���
	void print_obstacle(int y, int x){
		for(int i = 0; i < data[y][x].size(); i++){
			cout << "obstacle: "<< data[y][x][i].x << " "<<data[y][x][i].y<< endl;
		}
	}

	//������
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
			//�s��̗v�f�����v�Z
			int xnum = (int)((xmax - xmin) / rmax) + 1;
			int ynum = (int)((ymax - ymin) / rmax) + 1;
			//�z��̃T�C�Y�ύX
			data.resize(ynum);		// ()���̐������v�f���ɂȂ�
			for(int i = 0; i < ynum; i++ ){
				data[i].resize(xnum);
			}
			//�f�[�^�ǂݍ���
			std::ifstream ifs(str);
			if (!ifs) {
				cerr << "file failure" << endl;
			}
			while (!ifs.eof()){
				double xin,yin;
				//�Ƃ肠��������
				ifs >> xin >> yin;
				//�ǂ̋��ɓ��邩�v�Z
				int y = (int)((yin - ymin) / rmax);
				int x = (int)((xin - xmin) / rmax);
				//���Y����push
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

//erase���g�����߂̏���
template<typename T>
void remove(std::vector<T>& vector, unsigned int index){
	vector.erase(vector.begin() + index);
}
//�ǂɋ߂����闣�U�_������
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

// �z����ɓ������W�����݂���ꍇ��push_back���Ȃ����߂̊֐�
// �������W�ɗ��U�_����낤�Ƃ��Ă����ꍇ�͔��a���傫���ق����폜
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

//��Q���ߖT�̑傫���m�[�h�𗣎U��
void make_small_node(MAP *map, double s_x, double s_y, double s_r, double wall_x, double wall_y){
  double s_xx, s_yy, s_rr;
  double r = log2(map->dest[0].dr / s_r);
  const double dx[5] = {0, -1, +1, -1, +1};
  const double dy[5] = {0, +1, +1, -1, -1};
  const double dz[4][2] = {{+1, 0}, {-1, 0}, {0, +1}, {0, -1}};
  double dz_x, dz_y;

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
        && point_len(dz_x, wall_x, dz_y, wall_y) <  s_r * s_r
        && point_len(dz_x, wall_x, dz_y, wall_y) > s_rr * s_rr
        && r < density - 1){
					if(check_same(map, dz_x, dz_y, s_rr) == 0)
					map->dest.push_back(DEST(dz_x, dz_y, s_rr));
				}
    }
		if(s_rr < map->dest[0].dr * pow(0.5, density)) break;//�\���ȍׂ���������������I��
    else if(point_len(s_xx, wall_x, s_yy, wall_y) < s_rr * s_rr)//��Q�����m�[�h���Ɋ܂ނƂ�,���ׂ������U���D
      make_small_node(map, s_xx, s_yy, s_rr, wall_x, wall_y);
    else{//��Q�����m�[�h���Ɋ܂܂Ȃ�����,���U�_�Ƃ���D
			if(check_same(map, s_xx, s_yy, s_rr) == 0)
			map->dest.push_back(DEST(s_xx, s_yy, s_rr));
		}
  }
}

//�ǂƋ߂��Ƃ�1��Ԃ�
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

//���ʑS�̗̂��U��
void desti_all2(MAP *map, OBSTACLES obstacle, double xin, double yin){
	if(check_wall(map, obstacle, xin, yin, rin, 0) == 0){
		if(check_same(map, xin, yin, rin) == 0)
		map->dest.push_back(DEST(xin, yin, rin));
		}
		else map->stak.push_back(STAK(xin, yin, rin));
}

//���ʑS�̗̂��U��
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

//���U�����Ǘ�����֐�
void discretization(MAP *map, OBSTACLES obstacle){
	desti_all(map, obstacle);
	for(int s = 0; s < map->stak.size(); s++){
		check_wall(map, obstacle, map->stak[s].sx, map->stak[s].sy, map->stak[s].sr, 1);
	}
	erase_wall(map, obstacle);
}
//�f�[�^�C���|�[�g
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
	//�����l
  //cin >>x >>y;
	//cin >>n;
	//�Ǎ��W�C���|�[�g
	OBSTACLES obstacle("test_dots.txt", xmin, xmax, ymin, ymax, rmax);
	//obstacle.print_obstacle(5, 5);//�Ώۂ̃O���b�h�ɓ����Q���ꗗ�\��
	//cin >>gx >>gy;
	timer.start();
	discretization(&map, obstacle);
	timer.end_print();
	map.print();
}
