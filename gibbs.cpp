/* 
 * Gibbs sampler
 *
 * Student Name: Anbang Hu
 * Student AndrewID: anbangh
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <stdlib.h>
#include <cstdio>
#include <string.h>
#include <math.h>
#include <ctime>

using namespace std;

/* Function declarations */
void load_data(int **&, int &);
void initialization(int ***&, double **&, double **&, int, int);
void count_links(int **&, int ***&, int **&, int ***&, int, int);
void sample_target(int ***&, int **&, int **&, int ***&, int, int, double, double *);
void compute_posteriors(double **&, double **&, int **&, int ***&, int, int, double, double *);
double compute_objective(int **&, double **&, double **&, int, int);
void display(double **, int);
void save_objectives(double *, int);
void new_data1d(int *&y, int);
void new_data1d(double *&, int);
void new_data2d(int **&, int, int);
void new_data2d(double **&, int, int);
void new_data3d(int ***&, int, int, int);
void new_data3d(double ***&, int, int, int);
void free_data1d(int *&y, int);
void free_data1d(double *&, int);
void free_data2d(int **&, int, int);
void free_data2d(double **&, int, int);
void free_data3d(int ***&, int, int, int);
void free_data3d(double ***&, int, int, int);

// Driver
int main() {
	// Load data from hw3train.data
	int N;
	int **y;
	load_data(y, N);
	
	// Initialization
	int K 				= 5;
	double alpha 		= 0.02;
	double eta[] 		= {0.01, 0.05};
	int T 				= 10000;
	double obj[10000] 	= {0};

	int ***z;
	double **theta, **beta;
	initialization(z, theta, beta, N, K);

	// Display parameter information
	cout << "K     = " << K << endl;
	cout << "alpha = " << alpha << endl;
	cout << "eta0  = " << eta[0] << endl;
	cout << "eta1  = " << eta[1] << endl;
	cout << "T     = " << T << endl;

	// Count links
	int **node_comm, ***comm_comm;
	count_links(node_comm, comm_comm, y, z, N, K);

	// Sampling
	clock_t start = clock();
	for (int t = 0; t < T; t++) {
		// Sample target latent variables z
		sample_target(z, y, node_comm, comm_comm, N, K, alpha, eta);

		// Compute posterior theta and beta
		compute_posteriors(theta, beta, node_comm, comm_comm, N, K, alpha, eta);

		// Compute objective
		obj[t] = compute_objective(y, theta, beta, N, K);

		// Print epoch information
		if (t % 100 == 0) {
			cout << "[Iteration " << t << "]:\t Objective = " << obj[t];	
			cout << "\t(Time elapsed: " << (clock() - start ) / (double)CLOCKS_PER_SEC << " seconds)" << endl;
		}
	}

	// Save objectives
	save_objectives(obj, N);

	// Print beta
	display(beta, K);

	// Free memory
	free_data2d(y, N, N); 
	free_data3d(z, N, N, 2);
	free_data2d(theta, N, K);
	free_data2d(beta, K, K);
	free_data2d(node_comm, N, K);
	free_data3d(comm_comm, K, K, 2);

	return 0;
}

/* Function definitions */
void load_data(int **&y, int &N) {
	// Open file
	ifstream infile("./hw3train.data");
	if (!infile.is_open()) {
		cerr << "./hw3train.data cannot be opened" << endl;
	}
	// Obtain number of vertices
	string line;
	getline(infile, line);
	N = stoi(line);
	new_data2d(y, N, N);
	for (int i = 0; i < N; i++) {
		getline(infile, line);
		istringstream iss(line);
		vector<string> tokens;
		copy(istream_iterator<string>(iss),
		     istream_iterator<string>(),
		     back_inserter(tokens));
        for (unsigned int j = 0; j < tokens.size(); j++) {
        	y[i][stoi(tokens[j])] = 1;
        }
	}
	// Close file
	infile.close();
}

void initialization(int ***&z, double **&theta, double **&beta, int N, int K) {
	// Initialize z
	new_data3d(z, N, N, 2);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int w = 0; w < 2; w++) {
				z[i][j][w] = rand() % K;
			}
		}
	}

	// Initialize theta
	new_data2d(theta, N, K);
	double *tmp_theta_sum;
	new_data1d(tmp_theta_sum, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < K; j++) {
			theta[i][j] = 1.f * rand();
			// theta[i][j] = 1.f * (rand() % (N * K));
			tmp_theta_sum[i] += theta[i][j];
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < K; j++)
			theta[i][j] /= tmp_theta_sum[i];
	free_data1d(tmp_theta_sum, N);

	// Initialize beta
	new_data2d(beta, K, K);
	for (int i = 0; i < K; i++)
		for (int j = 0; j < K; j++)
			beta[i][j] = 1.f * rand() / RAND_MAX;
}

void count_links(int **&node_comm, int ***&comm_comm, int **&y, int ***&z, int N, int K) {
	// Count the number of links for node-community
	new_data2d(node_comm, N, K);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) continue;
			// Increment corresponding bucket
			node_comm[i][z[i][j][0]]++;
			node_comm[i][z[i][j][1]]++;
		}
	}

	new_data3d(comm_comm, K, K, 2);
	// Count the number of links for negative/positive community-community
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) continue;
			int direction;
			if (i < j) { // outgoing direction
				direction = 0;
			} else { // incoming direction
				direction = 1;
			}
			comm_comm[z[i][j][direction]][z[j][i][direction]][y[i][j]]++;
		}
	}
}

void sample_target(int ***&z, int **&y, int **&node_comm, int ***&comm_comm, int N, int K, double alpha, double *eta) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) continue;

			// Take out links that should not be included
			int direction = i < j ? 0 : 1;
			node_comm[i][z[i][j][direction]]--;
			node_comm[j][z[j][i][direction]]--;
			comm_comm[z[i][j][direction]][z[j][i][direction]][y[i][j]]--;

			// Compute prob distribution from which to sample
			double **p;
			new_data2d(p, K, K);
			double sum = 0;
			for (int k_i = 0; k_i < K; k_i++) {
				for (int k_j = 0; k_j < K; k_j++) {
					p[k_i][k_j] = (node_comm[i][k_i]+alpha)*(node_comm[j][k_j]+alpha)
									*(comm_comm[k_i][k_j][y[i][j]]+eta[y[i][j]])
									/(comm_comm[k_i][k_j][0] + comm_comm[k_i][k_j][1] + eta[0] + eta[1]);
					sum += p[k_i][k_j];
				}
 			}
 			for (int k_i = 0; k_i < K; k_i++) {
 				for (int k_j = 0; k_j < K; k_j++) {
 					p[k_i][k_j] /= sum;
 				}
 			}
 			sum = 0;
 			for (int w = 0; w < K*K; w++) {
 				int k_i = w / K;
 				int k_j = w % K;
 				sum += p[k_i][k_j];
 				p[k_i][k_j] = sum;
 			}

 			// Sample point
 			double r = 1.f * rand()/RAND_MAX;
 			bool success = false;
 			int z1, z2;
 			for (int w = 0; w < K*K; w++) {
 				int k_i = w / K;
 				int k_j = w % K;
 				if (p[k_i][k_j] > r) {
 					success = true;
 					z1 = k_i;
 					z2 = k_j;
 					break;
 				}
 			}
 			if (!success) {
 				z1 = K-1;
 				z2 = K-1;
 			}

 			// Update variables
 			node_comm[i][z1]++;
 			node_comm[j][z2]++;
 			comm_comm[z1][z2][y[i][j]]++;
 			z[i][j][direction] = z1;
 			z[j][i][direction] = z2;

 			// Free temporary memory
 			free_data2d(p, K, K);
		}
	}
}

void compute_posteriors(double **&theta, double **&beta, int **&node_comm, int ***&comm_comm, int N, int K, double alpha, double *eta) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < K; j++) {
			theta[i][j] = (node_comm[i][j] + alpha)/(2*(N-1) + K*alpha);
		}
 	}

 	for (int i = 0; i < K; i++) {
 		for (int j = 0; j < K; j++) {
 			beta[i][j] = (comm_comm[i][j][1] + eta[1])/(comm_comm[i][j][1] + comm_comm[i][j][0] + eta[1] + eta[0]);
 		}
 	}
}

double compute_objective(int **&y, double **&theta, double **&beta, int N, int K) {
	double lld = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) continue;
			double curr_lld = 0;
			for (int k_i = 0; k_i < K; k_i++) {
				for (int k_j = 0; k_j < K; k_j++) {
					if (y[i][j]) curr_lld += theta[i][k_i]*theta[j][k_j]*beta[k_i][k_j];
					else 	     curr_lld += theta[i][k_i]*theta[j][k_j]*(1-beta[k_i][k_j]);
				}
			}
			lld += log(curr_lld);
		}
	}
	return lld;
}

void save_objectives(double *obj, int N) {
	ofstream myfile ("objective.txt");
	for (int i = 0; i < N; i++)
		myfile << obj[i] << endl;
}

void display(double **y, int N) {
	cout << "y:" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << y[i][j] << ",";
		}
		cout << endl;
	}
}

void new_data1d(int *&y, int m) {
	y = new int[m];
	memset(y, 0, sizeof(int)*m);
}

void new_data2d(int **&y, int m, int n) {
	y = new int*[m];
	for (int i = 0; i < m; i++) {
		y[i] = new int[n];
		memset(y[i], 0, sizeof(int)*n);
	}
}

void new_data3d(int ***&y, int m, int n, int k) {
	y = new int**[m];
	for (int i = 0; i < m; i++)
		y[i] = new int*[n];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			y[i][j] = new int[k];
			memset(y[i][j], 0, sizeof(int)*k);
		}
	}
}

void new_data1d(double *&y, int m) {
	y = new double[m];
	memset(y, 0, sizeof(double)*m);
}

void new_data2d(double **&y, int m, int n) {
	y = new double*[m];
	for (int i = 0; i < m; i++) {
		y[i] = new double[n];
		memset(y[i], 0, sizeof(double)*n);
	}
}

void new_data3d(double ***&y, int m, int n, int k) {
	y = new double**[m];
	for (int i = 0; i < m; i++)
		y[i] = new double*[n];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			y[i][j] = new double[k];
			memset(y[i][j], 0, sizeof(double)*k);
		}
	}
}

void free_data1d(int *&y, int m) {
	delete [] y;
}

void free_data2d(int **&y, int m, int n) {
	for (int i = 0; i < m; i++) {
		delete [] y[i];
	}
}

void free_data3d(int ***&y, int m, int n, int k) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			delete [] y[i][j];
		}
	}
	for (int i = 0; i < m; i++) {
		delete [] y[i];
	}
}

void free_data1d(double *&y, int m) {
	delete [] y;
}

void free_data2d(double **&y, int m, int n) {
	for (int i = 0; i < m; i++) {
		delete [] y[i];
	}
}

void free_data3d(double ***&y, int m, int n, int k) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			delete [] y[i][j];
		}
	}
	for (int i = 0; i < m; i++) {
		delete [] y[i];
	}
}