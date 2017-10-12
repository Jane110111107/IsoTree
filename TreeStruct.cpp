#include "TreeStruct.h"
#include<iostream>
#include<vector>

using namespace std;


TreeNode::TreeNode(int id, float cov) {
	
	treenode_id = id;
	treenode_cov = cov;

}


bool TreeNode::add_child(int child) {

	if (child < 0)
		return false;

	if(!children.empty()) {
		
		for (size_t i = 0; i < children.size(); i++) 
			if (children[i] == child)
				return false;
	}

	this -> children.push_back(child);
	return true;

}


void TreeNode::add_candi(int p, float x){

	bool is_add = true;
	int add_id = -1;
	for (int i = 0; i < candi.size(); i++) {
		if (candi[i] == p){
			is_add = false;
			add_id = i;
			break;
		}
	}
	if (is_add) {
		candi.push_back(p);
		priority.push_back(x);
	} else {
		priority[add_id] = priority[add_id] + x;
	}
}


void TreeNode::dele_candi(int p){

	for (int i = 0; i < candi.size(); i++) {
		if (candi[i] == p) {
			candi.erase(candi.begin()+i);
			priority.erase(priority.begin()+i);
			break;
		}
	}

}


TreeStruct::TreeStruct() {

	treenode_sum = 0;

}


int TreeStruct::add_node(TreeNode node) {

	treenode_set.push_back(node);
	treenode_sum++;
	return treenode_sum-1;

}

int TreeStruct::total_node() {

	return treenode_sum;

}


void TreeStruct::add_child(int v, int child) {

	treenode_set[v].add_child(child);

}


void TreeStruct::dfs_trans(vector<float>& w_trans, vector<vector<int> >& trans_vec) {

	w_trans.clear();
	trans_vec.clear();
	vector<int> temp;

	dfs_v(0,w_trans,temp, trans_vec);

}


void TreeStruct::dfs_v(int node_id, vector<float>& w_trans, vector<int>& temp, vector<vector<int> >& trans_vec) {

	int child_sum = treenode_set[node_id].children.size();
	//cout << child_sum << endl;////
	temp.push_back(treenode_set[node_id].treenode_id);
	//cout << temp.size() << endl;////
	if (child_sum == 0) {

		trans_vec.push_back(temp);
		w_trans.push_back(treenode_set[node_id].treenode_cov);
		temp.pop_back();

	} else {

		for (int i = 0; i < child_sum; i++)
			dfs_v(treenode_set[node_id].children[i], w_trans, temp, trans_vec);

		temp.pop_back();

	}

}
