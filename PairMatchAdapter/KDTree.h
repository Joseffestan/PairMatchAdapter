#pragma once
#include "pfh.h"
/* the maximum number of keypoint NN candidates to check during BBF search */
#define KDTREE_BBF_MAX_NN_CHKS 200  

/* threshold on squared ratio of distances between NN and 2nd NN */
#define NN_SQ_DIST_RATIO_THR 0.49 

//kd tree节点PFH
struct KDNodePoint {
	int ki;                      /**< partition key index */
	double kv;                   /**< partition key value */
	int leaf;                    /**< 1 if node is a leaf, 0 otherwise */
	PFH* pfh;    /**< points at this node */
	int n;                       /**< number of points */
	struct KDNodePoint* kd_left;     /**< left child */
	struct KDNodePoint* kd_right;    /**< right child */
};
struct bbf_data
{
	double d;
	void* old_data;
};
//创建KD Tree 
extern struct KDNodePoint* KDTreeBuild(PFH* points, int n);
//KD Tree搜索最邻近 
extern int KDTreeSearch(KDNodePoint* kd_root, PFH* point, int k, PFH*** nbrs, int max_nn_chks);
//KD Tree释放
extern void KDTreeRelease(KDNodePoint* kd_root);