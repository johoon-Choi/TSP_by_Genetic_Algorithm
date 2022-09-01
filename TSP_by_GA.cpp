// Vertex 0 1 2 3 4 ..
//        A B C D E ..\\


// ���յ��� �����ϴ� ������ �߿��ϴ�. ������.

#include <iostream>
#include <iomanip>
#include <queue>
#include <string>
#include <cstdlib>
#include <ctime>
#include <set>

#define SEED_KEY 1000
#define INF 1000000
#define N 26
#define POPULATION_NUM 300
#define START_VERTEX 0
#define DISTANCE_RANGE 50
#define TEST_CASE 100
#define K 5				// ���յ� ���Ŀ� �ʿ�
						// ���� 1�� ����� �� ���� ���յ����� Ŀ����. ���� �о�����.
#define Pc 700			// Cross Over �� Ȯ�� : 70%�� �ǹ�
#define Pm 1			// Mutation �� Ȯ�� : 0.1%�� �ǹ�

template <typename E>
inline void Swap(E& a, E& b) {
	int temp = a;
	a = b;
	b = temp;
}

inline int fit_func(double Fmax, double Fmin, double Fcur) {
	return (Fmax - Fcur) + (Fmax - Fmin) / (K - 1.0);
} // http://www.aistudy.com/biology/genetic/operator_moon.htm ����

using namespace std;

int max_W = 0;							// Selection �۾��� �����ȴ�.
int min_W = INF;
int genes_fit_sum = 0;
const int start_vertex = START_VERTEX;

struct population {
	string gene = "";
	int fit = 0;
};

void Initialization(int** W, population* parents, const int n) {
	string cur_gene = "";
	int visited = 0;
	int before = ::start_vertex;
	int cur_dis_sum  = 0;
	int max_dis_sum = 0;
	int min_dis_sum = INF;

	srand((unsigned)time(NULL));

	for (int i = 0; i < POPULATION_NUM; i++) {
		while (visited != (1 << n) - (1 + (1 << start_vertex))) {
			// start_vertex == 3 , n == 5 ?
			// 100000(2) - 1 - (1 << 3) == 10111(2) ���� ����
			int cur = START_VERTEX;		// [0, n - 1] ������ ���� ����
			while (cur == START_VERTEX) cur = rand() % n;

			if (visited & (1 << cur)) continue;
			
			cur_dis_sum += W[before][cur];			// ���յ� �κп� before --- cur �� �Ÿ����� �߰�.

			cur_gene += char(cur + 'A');		// cur_gene �� ���� vertex�� �ش�Ǵ� ���⼭�� �߰�
			before = cur;						// ���� vertex üũ ���� ���� vertex�� ���� vertex�� ���
			visited += (1 << cur);				// ���� vertex �湮 ó��
		}

		cur_dis_sum += W[before][start_vertex];

		parents[i].gene = cur_gene;
		parents[i].fit = cur_dis_sum;

		max_dis_sum = max(max_dis_sum, cur_dis_sum);
		min_dis_sum = min(min_dis_sum, cur_dis_sum);

		cur_gene = "";		// ���� �۾��� ���� gene_temp�� �ٽ� �ʱ�ȭ
		visited = 0;		// ���� �۾��� ���� visited�� �ٽ� �ʱ�ȭ
		before = 0;			// ���� �۾��� ���� before�� �ٽ� �ʱ�ȭ
		cur_dis_sum = 0;		// ���� �۾��� ���� edge_sum�� �ٽ� �ʱ�ȭ
	}

	//cout << " max fit sum is " << max_dis_sum << '\n';

	for (int i = 0; i < POPULATION_NUM; i++) {
		int cur_dis_sum = parents[i].fit;
		parents[i].fit = fit_func(max_dis_sum, min_dis_sum, cur_dis_sum);
		::genes_fit_sum += parents[i].fit;
	}
}

void Selection(population* parents, const int n, int* S) {
	int rand_fit[2] = {}; 

	for (int i = 0; i < 2; i++) {
		rand_fit[i] = rand() % (::genes_fit_sum + 1);

		int fit_sum_temp = 0; // �˸��� ��ġ�� gene�� get �ϱ� ����
		int j = 0;

		while (1) {
			fit_sum_temp += parents[j].fit;
			if (rand_fit[i] < fit_sum_temp) break;
			j++;
		}

		S[i] = j; // parents ���� j ��° gene �� index�� j �� S[i] �� ����
	}
}

void CrossOver(int** W, int n, population* parents, int* S) {
	int s1 = S[0];
	int s2 = S[1];

	population* p1 = &parents[s1];
	population* p2 = &parents[s2];

	const int len = N;
	int left = rand() % len;
	int right = left;
	while (left == right) right = rand() % len;

	if (right < left) Swap(left, right);

	set<char> set1, set2;
	// ############# �������� Cross Over ���� #############
	for (int i = left; i <= right; i++) {
		Swap((*p1).gene[i], (*p2).gene[i]);
		set1.insert((*p1).gene[i]);
		set2.insert((*p2).gene[i]);
	}

	string g1 = (*p1).gene;
	string g2 = (*p2).gene;
	 
	queue<int> q1, q2;	// ���� gene �� ���� ���⼭���� �ִ´�.

	for (int i = 0; i < N; i++) {
		if (!set1.count(g1[i])) set1.insert(g1[i]);
		else {
			q1.push(g1[i]);
		}

		if (!set2.count(g2[i])) set2.insert(g2[i]);
		else {
			q2.push(g2[i]);
		}
	}
}

int TSP(const int n, int** W) {
	int selected[2];									// ���õ� gene1, gene2 �� index value
	population* parents = new population[POPULATION_NUM];

	Initialization(W, parents, n);				// ������ �����ڵ� ���� �Ϸ�

	for (int i = 0; i < POPULATION_NUM; i++) {
		cout << "gene[" << i << "] = " << parents[i].gene << " : " << parents[i].fit << '\n';
	}

	int tc = 0;

	int c_cnt = 0;
	int m_cnt = 0;

	while (tc++ < TEST_CASE) {
		// ############## SELECTION ##############
		Selection(parents, n, selected);			// s1, s2

		// ############## CROSS OVER ##############
		if ((rand() % 1000) + 1 <= Pc) c_cnt++; // CrossOver Function

		// ############## MUTATION ##############
		if ((rand() % 1000) + 1 <= Pm) m_cnt++; // Mutation Function
	}

	cout << " Crossover count is " << c_cnt << '\n';
	cout << " Mutation count is " << m_cnt << '\n';

	return 0;
}

int main() {
	srand(SEED_KEY);

	int n = N;
	int** W = new int* [n];
	for (int i = 0; i < n; i++)	W[i] = new int[n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i != j) {
				W[i][j] = rand() % DISTANCE_RANGE + 1;
				if (min_W > W[i][j]) min_W = W[i][j];
				if (W[i][j] > max_W && W[i][j] != INF) max_W = W[i][j];
			}
			else W[i][j] = 0;

	clock_t start, end;
	start = clock();
	int ans = TSP(n, W);
	end = clock();
	
	cout << " The answer is " << ans << '\n' << '\n';
	cout << " �ҿ� �ð� : " << end - start << " ms" << '\n';

	for (int i = 0; i < n; i++) delete[] W[i];
	delete[] W;
}