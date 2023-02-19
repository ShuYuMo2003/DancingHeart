#include <graphics.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <Windows.h>
#include <vector>
#include <cstdlib>
#include <thread>

#pragma comment(lib,"winmm.lib")
#define debug() (std::cerr << "Running at " << __LINE__ << std::endl)
using namespace std;

#define BUFFER_SIZE 130
int MainFrameId = 0;

namespace Beating_Heart {
	using std::cerr;
	using std::endl;
	using std::vector;

	typedef double db;

	const int WIDTH = 1280;
	const int HEIGH = 720;
	const int CNSH = 30;
	const db PI = acos(-1.0);
	const db sigma = 0.6;
	const int RNAD_SHIFT = 30;


	const COLORREF MAIN_COLOR = RGB(0xFF, 0x9E, 0x9E);

	const int W_L = -(WIDTH >> 1);
	const int W_U = (WIDTH >> 1);
	const int H_L = -(HEIGH >> 1);
	const int H_U = (HEIGH >> 1);

	COLORREF _graph[WIDTH][HEIGH];

	struct Point {
		int x, y;
		Point(int _x, int _y) { x = _x; y = _y; }
		Point() {};
		Point operator - (const Point& rhs) { return Point(x - rhs.x, y - rhs.y); }
		Point operator + (const Point& rhs) { return Point(x + rhs.x, y + rhs.y); }
		Point operator * (const db& t) { return Point(x * t, y * t); }
		double len() { return sqrt((db)x * x + (db)y * y); }
		double slope() { return (double)y / x; }
		void print() { printf("(%d, %d)\n", x, y); }
	};

	struct Pixel {
		Point pos;
		COLORREF c;
		Pixel(const Point& $0, const COLORREF& $1) { pos = $0; c = $1; }
	};

	template<typename T> T _cube(T x) { return x * x * x; }
	template<typename T> T _square(T x) { return x * x; }
	template<typename T> T _max(T a, T b) { return a > b ? a : b; }
	template<typename T> T _min(T a, T b) { return a < b ? a : b; }
	template<typename T> T _abs(T a) { return a < 0 ? -a : a; }
	template<typename T> T _lerp(T a, T b, T radio) { return (a + (b - a) * radio); }



	#define GUASS_SIZE 5000
	db gauss_value[GUASS_SIZE];

	db gauss_function_(int len_2) {
		return (1. / (2. * PI * sigma * sigma)) * exp(-len_2 / (2 * sigma * sigma));
	}

	void initialization() {
		initgraph(WIDTH, HEIGH);
		for (int i = 0; i < GUASS_SIZE; i++) gauss_value[i] = gauss_function_(i);
	}

	COLORREF temp;
	COLORREF& _(Point now) {
		int x = now.x + (WIDTH >> 1), y = -now.y + (HEIGH >> 1) + CNSH;
		if (x < 0 || y < 0 || x > WIDTH || y > HEIGH) return temp;
		return _graph[x][y];
	}

	bool _rand(double p) { // 返回 true 的概率是 p.
		assert(p <= 1 && p >= 0);
		return rand() % int(1 / p) == 0;
	}

	double randdb(db L, db R) { // 随机返回 0-1
		return (rand() / (double)RAND_MAX) * (R - L) + L;
	}

	db heart_function(Point now, db sc = 220) {
		db x = now.x / sc;
		db y = now.y / sc;
		return _cube(x * x + y * y - 1) - x * x * y * y * y;
	}

	COLORREF genColor(db scales) {
		const db hue = 0.999;
		COLORREF ret = (scales >= 1 ? HSLtoRGB(hue, 2 - scales, 1) : HSLtoRGB(hue, 1, scales));
		return ret;
		// return RGB(_min(GetRValue(ret) / 256. * 1.5, 255.),
		//            _min(GetGValue(ret) / 256. * 1.5, 255.),
		//            _min(GetBValue(ret) / 256. * 1.5, 255.));
	}

	COLORREF muticolor(COLORREF now, db scale) {
		return RGB((int)(GetRValue(now) * scale),
			(int)(GetGValue(now) * scale),
			(int)(GetBValue(now) * scale));
	}

	void flushGraph() {
		static COLORREF old[WIDTH][HEIGH] = { 0 };
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < HEIGH; j++) {
				if (old[i][j] != _graph[i][j]) {
					old[i][j] = _graph[i][j];
					putpixel(i, j, old[i][j] = _graph[i][j]);
				}
			}
		}
	}

	void clearPraph() { memset(_graph, 0, sizeof(_graph)); }

	void gaussFilter(COLORREF _graph[WIDTH][HEIGH]) {
		static COLORREF src[WIDTH][HEIGH];
		memset(src, 0, sizeof(src));

		for (int i = 0; i < WIDTH; i++)
			for (int j = 0; j < HEIGH; j++) {
				if (!_graph[i][j]) // 跳过纯黑像素.
					continue;
				Point _0(_max(0., i - 3 * sigma), _max(0., j - 3 * sigma));
				Point _1(_min(WIDTH - 1., i + 3 * sigma), _min(HEIGH - 1., j + 3 * sigma));
				for (int x = _0.x; x <= _1.x; x++)
					for (int y = _0.y; y <= _1.y; y++) {
						src[x][y] += muticolor(_graph[i][j], gauss_value[_square(x - i) + _square(y - j)]);
					}
			}
		memcpy(_graph, src, sizeof(src));
	}

	void render(const vector<Pixel>& P) {
		for (const Pixel& now : P) {
			// Point _0(now.pos.x - 3 * sigma, now.pos.y - 3 * sigma);
			// Point _1(now.pos.x + 3 * sigma, now.pos.y + 3 * sigma);
			_(now.pos) = now.c;
			// for(int i = _0.x; i <= _1.x; i++) {
			//     for(int j = _0.y; j <= _1.y; j++) {
			//         _(Point(i, j)) += muticolor(now.c, gauss_value[_square(i - now.pos.x) + _square(j - now.pos.y)]);
			//     }
			// }
		}
	}

	vector<Pixel> generateSimpleInnerHeart(int pointsCnt) {
		Point center(0, CNSH);
		static Point POINTS_POOL[WIDTH * HEIGH];
		vector<Pixel> ret;
		ret.clear();
		int POOL_LEN = 0;

		for (int i = W_L; i < W_U; i++) {
			for (int j = H_L; j < H_U; j++) { // 生成边界点
				if (heart_function(Point(i, j), 220) * heart_function(Point(i, j), 216) < 0) {
					POINTS_POOL[POOL_LEN++] = Point(i, j);
				}
			}
		}
		printf("Unique points total = %d\n", POOL_LEN);
		int repeat = _max(1, pointsCnt / POOL_LEN);
		for (int i = 0; i < POOL_LEN; i++) { // 向内扩散
			Point now = POINTS_POOL[i];
			for (int t = 0; t < repeat; t++) {
				db shift = _max(1 - 0.7 * 0.25 * (-log(randdb(0, 1))), 0.);
				ret.push_back(Pixel((now - center) * shift + center,
					genColor(pow(shift * 0.8 * randdb(0.65, 1.2), 1))));
			}
		}
		for (auto& now : ret) {// 扣去中心异常密集点
			db d = fabs((now.pos - center).len());
			if ((d < 4) || (d < 9 && _rand(0.5)))
				now.c = 0;
		}
		return ret;
	}

	vector<Pixel> generateSimpleOuterHeart(db rho) {
		vector<Pixel> ret;
		for (int i = W_L; i < W_U; i++) {
			for (int j = H_L; j < H_U; j++) {
				db shift = _max(1 - 0.7 * 0.25 * (-log(randdb(0, 1))), 0.);
				if (heart_function(Point(i, j), 290) < 0 && _rand(rho))
					ret.push_back(Pixel(Point(i, j),
						genColor(shift * 0.5 * randdb(0.9, 1.2))));
			}
		}
		return ret;
	}

	void randomShift(vector<Pixel>& P, int RNAD_SHIFT) {
		for (auto& now : P) {
			now.pos.x += randdb(-RNAD_SHIFT, RNAD_SHIFT);
			now.pos.y += randdb(-RNAD_SHIFT, RNAD_SHIFT);
		}
	}


	vector<Pixel> sphereWithRandomShift(vector<Pixel>& P, db radio, int RAND_SHIFT = 1) {
		Point center(0, CNSH);
		vector<Pixel> ret; ret.clear();
		db radius = 10;
		for (auto now : P) {
			now.pos = now.pos * _lerp(1., 1.5, radio);

			db dist = (now.pos - center).len();

			now.pos = now.pos * _lerp(0.9, radius / dist, radio * 0.4);

			now.pos.x += (rand() % 2 - 1);
			now.pos.y += (rand() % 2 - 1);
			ret.push_back(now);
		}

		return ret;
	}

	IMAGE buffer[BUFFER_SIZE];

	void PreProcessHeart() {
		srand((unsigned)time(0));
		// initialization();

		vector<Pixel> _inner = generateSimpleInnerHeart(90000);
		vector<Pixel> _outer = generateSimpleOuterHeart(0.04);
		randomShift(_outer, RNAD_SHIFT);

		db Rad_L = 0.8;
		db Rad_R = 1.2;

		int L = 0, R = BUFFER_SIZE - 1;
		double flame = 2. / BUFFER_SIZE, rad = 0;

		for (int t = 0; t < BUFFER_SIZE; t++) {
			db radio = _square(rad) * ((Rad_R - Rad_L) / 4) + Rad_L;
			rad += flame;
			clearPraph();
			render(sphereWithRandomShift(_inner, radio));
			render(sphereWithRandomShift(_outer, 0.8 * ((Rad_R - radio) + Rad_L)));
			gaussFilter(_graph);
			flushGraph();
			getimage(&buffer[t], 0, 0, WIDTH, HEIGH);
			printf("\rProcessing frame:"); int __end = 100 * double(t + 1) / BUFFER_SIZE;
			for (int i = 0; i < __end; i++) putchar('=');
		}
		printf("\n");
	}

	void flush(int id) {
		id = max(0, id);
		id = min(BUFFER_SIZE, id);
		static int oldid = 0;
		if (oldid == id) return;
		putimage(0, 0, &buffer[oldid = id]);
		printf("\r\rscale:", MainFrameId);
		for (int i = 0; i < (BUFFER_SIZE - MainFrameId) / 1.5; i++) putchar('=');
		for (int i = 0; i < MainFrameId / 1.5; i++) putchar(' ');
		// cerr << id << endl;
	}
	void handle() { while (true) flush(MainFrameId); }
}


const int VOLUMN = 100;
const int FPM = 70;
const double FULLNOTE = 4 * 60 * (double)1000 / FPM;
HMIDIOUT MainHandle;
const int FPS = 40;
#define combine(a, b, c) (((int)(a) << 16) + ((int)(b) << 8) + (int)c)

enum Scale
{
	Rest = 0, C8 = 108,
	B7 = 107, A7s = 106, A7 = 105, G7s = 104, G7 = 103, F7s = 102, F7 = 101, E7 = 100, D7s = 99, D7 = 98, C7s = 97, C7 = 96,
	B6 = 95, A6s = 94, A6 = 93, G6s = 92, G6 = 91, F6s = 90, F6 = 89, E6 = 88, D6s = 87, D6 = 86, C6s = 85, C6 = 84,
	B5 = 83, A5s = 82, A5 = 81, G5s = 80, G5 = 79, F5s = 78, F5 = 77, E5 = 76, D5s = 75, D5 = 74, C5s = 73, C5 = 72,
	B4 = 71, A4s = 70, A4 = 69, G4s = 68, G4 = 67, F4s = 66, F4 = 65, E4 = 64, D4s = 63, D4 = 62, C4s = 61, C4 = 60,
	B3 = 59, A3s = 58, A3 = 57, G3s = 56, G3 = 55, F3s = 54, F3 = 53, E3 = 52, D3s = 51, D3 = 50, C3s = 49, C3 = 48,
	B2 = 47, A2s = 46, A2 = 45, G2s = 44, G2 = 43, F2s = 42, F2 = 41, E2 = 40, D2s = 39, D2 = 38, C2s = 37, C2 = 36,
	B1 = 35, A1s = 34, A1 = 33, G1s = 32, G1 = 31, F1s = 30, F1 = 29, E1 = 28, D1s = 27, D1 = 26, C1s = 25, C1 = 24,
	B0 = 23, A0s = 22, A0 = 21
};


struct Note_t {
	double absolute_time;
	int level_code;
	double length;
	int volumn;
	string desc;
	bool operator < (const Note_t& rhs) const { return absolute_time < rhs.absolute_time; }
	Note_t() {}
	Note_t(double ti) { absolute_time = ti; }
	Note_t(double a, int b, double c, int d, string e) {
		absolute_time = a;
		level_code = b;
		length = c;
		volumn = d;
		desc = e;
	}
};

int StringCnt(const string& src, char aim) {
	int ret = 0;
	for (auto now : src)
		ret += (aim == now);
	return ret;
}

vector <Note_t> noteDecoder(string mask) { // string mask: [音符对应数字] [+高音 / *低音 / 省略 中音] [-对应简谱延长] [_对应简谱下划线] [.附点]
	static int do_unit_shift = 12;
	static int idx_shift[] = { INT_MIN, 0, 2, 4, 5, 7, 9, 11 }; // idx[i]: 与 do (idx[1]) 的偏移量
	static char msg[200];

	mask += "!"; // 结束符

	vector<int> level;
	int handling = -1;
	for (auto _ : mask) {
		if (_ <= '9' && _ >= '0') {
			if (handling != -1)
				level.push_back(handling);
			handling = idx_shift[_ - '0'] + C4;
		}
		else if (_ == '+') handling += do_unit_shift;
		else if (_ == '*') handling -= do_unit_shift;
		else {
			assert(handling != -1);
			level.push_back(handling);
			break;
		}
	}
	int length = 4; // 4分音符
	switch (StringCnt(mask, '-')) {
	case 3:
		length = 1;
		break;
	case 1:
		length = 2;
		break;
	}

	switch (StringCnt(mask, '_')) {
	case 1:
		length = 8;
		break;
	case 2:
		length = 16;
		break;
	}
	double timel = FULLNOTE / (double)length;
	if (StringCnt(mask, '.'))
		timel *= 1.5;
	sprintf_s(msg, "Length = 1 /% 2d, mask = %s. %s", length, mask.c_str(), level.size() < 2 ? "" : "（和弦）");
	vector<Note_t> ret;
	for (int note : level) {
		ret.push_back(Note_t(-1, note, timel, 0, string(msg)));
	}
	return ret;
}

void play(Note_t note, int playertype) {
	static double lastPlayTime[16] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1};

	int layer = -1;
	for (int i = 0; i < 16; i++) {
		if (lastPlayTime[i] != note.absolute_time || lastPlayTime[i] < 0) {
			layer = i;
			break;
		}
	}

	assert(layer != -1);

	lastPlayTime[layer] = note.absolute_time;

	printf("\n\nPlaying note: % 3d, layer: % 2d, absolute time: %10.4lf.\n - Decoder Info: %s\n\n", note.level_code, layer, note.absolute_time, note.desc.c_str());
	midiOutShortMsg(MainHandle, 0x7BB0 + layer);
	midiOutShortMsg(MainHandle, combine(0, playertype, (0xC0 + layer)));
	if(note.level_code > 0)
		midiOutShortMsg(MainHandle, (note.volumn << 16) + ((note.level_code ) << 8) + 0x90 + layer);
}

clock_t nowtime(clock_t startTime) { return ((clock() - startTime) / (double)CLOCKS_PER_SEC) * 1000; }


int main() {
	Beating_Heart::initialization();
	system("pause");

	midiOutOpen(&MainHandle, 0, 0, 0, CALLBACK_NULL);
	FILE* NotesFile[2];
	if (fopen_s(&NotesFile[0], "notes_0.txt", "r") != 0 or fopen_s(&NotesFile[1], "notes_1.txt", "r") != 0)
		return printf("Error on opening notes file."), -1;

	for (int i = 0; i < 16; i++) {
		midiOutShortMsg(MainHandle, combine(0, 1, (0xC0 + i)));
	}

	vector<Note_t> seq[2], all;
	char tempmask[100];
	double timeStamp = 0;
	vector<double> flagstamps;

	for (int SEQNO = 0; SEQNO < 2; SEQNO++) {
		timeStamp = 0;
		int noTimeStamp = 0;
		while (fscanf_s(NotesFile[SEQNO], "%s", tempmask, 100) != EOF) {
			// cout << "read: " << tempmask << endl;
			if (tempmask[0] == '#') {
				if (SEQNO == 0)
					flagstamps.push_back(timeStamp);
				else {
					assert(fabs(timeStamp - flagstamps[noTimeStamp]) < 0.2);
					timeStamp = flagstamps[noTimeStamp++];
				}
			}
			else {
				vector<Note_t> slice = noteDecoder(tempmask);
				for (auto now : slice) {
					now.absolute_time = timeStamp;
					now.desc += SEQNO == 0 ? "Seq = 0" : "Seq = 1";
					now.volumn = SEQNO == 0 ? 110 : 80;
					seq[SEQNO].push_back(now);
					all.push_back(now);
				}
				timeStamp += slice[0].length;
			}
		}
	}
	sort(all.begin(), all.end());
	printf("Sequenclize notes done\n");

	vector<pair<double, double> > aniSc;
	double __step = (double)1 / FPS, __end = (all.end() - 1)->absolute_time + (all.end() - 1)->length;
	for (double now = 0; now <= __end; now += __step)
		aniSc.push_back(make_pair(now, 0));

	int totalNote = all.size();
	int noteCnt = 0;
	auto startptr = aniSc.begin();
	int limL = 5;
	for (auto slice : all) {
		double L = slice.absolute_time - limL;
		double R = slice.absolute_time + slice.length;
		double height = ((slice.volumn * 0.01) * (slice.level_code * 0.5));
		if (slice.level_code < 0) continue;
		printf("\rProccessing effect of note No.% 3d / %d (%s)", noteCnt++, totalNote, slice.desc.c_str());
		while (startptr < aniSc.end() && startptr->first < L)
			startptr++;
		auto aim = startptr;
		while (aim < aniSc.end() && aim->first <= R) {
			double distance = aim->first - slice.absolute_time;
			aim->second += height * (distance > 0 ? (1 - distance / slice.length) : (1 + distance / limL));
			aim++;
		}
	}


	printf("Animation(seq) preprocess done\n");
	Beating_Heart::PreProcessHeart();
	printf("Animation(frame) preprocess done\n");

	double maxAniRadio = 0, minAniRadio = INT_MAX;
	for (auto& slice : aniSc) {
		if (slice.second < 0) continue;
		maxAniRadio = max(maxAniRadio, slice.second);
		minAniRadio = min(minAniRadio, slice.second);
	}


	// thread(Beating_Heart::handle).detach();


	clock_t startTime = clock();
	#define NOW (nowtime(startTime))
	unsigned int ptr[2] = {0}, aniptr = 0;

	int skip = 0;
	while (ptr[0] < seq[0].size() || ptr[1] < seq[1].size()) {
		for (int SEQNO = 0; SEQNO < 2; SEQNO++) {
			while (ptr[SEQNO] < seq[SEQNO].size() && seq[SEQNO][ptr[SEQNO]].absolute_time <= NOW) {
				play(seq[SEQNO][ptr[SEQNO]], SEQNO == 0 ? 1 : 1);
				ptr[SEQNO] ++;
			}
		}

		while (aniptr < aniSc.size() && aniSc[aniptr].first <= NOW) {
			MainFrameId = BUFFER_SIZE - int(((aniSc[aniptr].second - minAniRadio) / (maxAniRadio - minAniRadio)) * BUFFER_SIZE);
			aniptr++;
		}
		Beating_Heart::flush(MainFrameId);
		// cerr << MainFrameId << endl;
		//Sleep(FULLNOTE / 64);
	}
	for (MainFrameId = 0; MainFrameId < BUFFER_SIZE; MainFrameId++) {
		Beating_Heart::flush(MainFrameId);
		Sleep(30);
	}

	return 0;
}