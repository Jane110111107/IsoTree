#ifndef TREESTRUCT_H
#define TREESTRUCT_H


#include<iostream>
#include<algorithm>
#include<vector>


using namespace std;


class TreeNode {

public:

	int treenode_id;
	float treenode_cov;
	vector<int> children;
	vector<int> candi;
	vector<float> priority;

	TreeNode(int id, float cov);
	bool add_child(int child);
	void add_candi(int p, float x);
	void dele_candi(int p);

};


class TreeStruct {

public:

	vector<TreeNode> treenode_set;
	int treenode_sum;
	TreeStruct();
	int add_node(TreeNode node);
	int total_node();
	void add_child(int v, int child);
	void dfs_trans(vector<float>& w_trans, vector<vector<int> >& trans_vec);
	void dfs_v(int node_id, vector<float>& w_trans, vector<int>& temp, vector<vector<int> >& trans_vec);

};



#endif
