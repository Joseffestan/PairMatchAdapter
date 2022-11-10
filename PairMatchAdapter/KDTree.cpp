
#include "KDTree.h"
#include "minpq.h"
#include<stdlib.h>

struct KDNodePoint* kd_node_init(struct PFH* points, int n);
void expand_kd_node_subtree(struct KDNodePoint* kd_node);
void assign_part_key(struct KDNodePoint* kd_node);
double median_select(double*, int);
double rank_select(double*, int, int);
void insertion_sort(double*, int);
int partition_array(double*, int, double);
void partition_nodes(struct KDNodePoint* kd_node);
struct KDNodePoint* explore_to_leaf(struct KDNodePoint* kd_node, struct PFH* point, struct min_pq* min_pq);
double descr_dist_sq(struct PFH* f1, struct PFH* f2);
int insert_into_nbr_array(struct PFH* point, struct PFH** nbrs, int n, int k);

struct KDNodePoint* KDTreeBuild(PFH* points, int n)
{
	struct KDNodePoint* kd_root;

	if (!points || n <= 0)
	{
		return NULL;
	}

	kd_root = kd_node_init(points, n);
	expand_kd_node_subtree(kd_root);

	return kd_root;
}
/*
Initializes a kd tree node with a set of features.  The node is not
expanded, and no ordering is imposed on the features.

@param features an array of image features
@param n number of features

@return Returns an unexpanded kd-tree node.
*/
struct KDNodePoint* kd_node_init(struct PFH* points, int n)
{
	struct KDNodePoint* kd_node;

	kd_node = (struct KDNodePoint*)malloc(sizeof(struct KDNodePoint));
	memset(kd_node, 0, sizeof(struct KDNodePoint));
	kd_node->ki = -1;
	kd_node->pfh = points;
	kd_node->n = n;
	return kd_node;
}
/*
Recursively expands a specified kd tree node into a tree whose leaves
contain one entry each.
@param kd_node an unexpanded node in a kd tree
*/
void expand_kd_node_subtree(struct KDNodePoint* kd_node)
{
	/* base case: leaf node */
	if (kd_node->n == 1 || kd_node->n == 0)
	{
		kd_node->leaf = 1;
		return;
	}

	assign_part_key(kd_node);
	partition_nodes(kd_node);

	if (kd_node->kd_left)
		expand_kd_node_subtree(kd_node->kd_left);
	if (kd_node->kd_right)
		expand_kd_node_subtree(kd_node->kd_right);
}
/*
Determines the descriptor index at which and the value with which to
partition a kd tree node's features.

@param kd_node a kd tree node
*/
void assign_part_key(struct KDNodePoint* kd_node)
{
	PFH* points;
	double kv, x, mean, var, var_max = 0;
	double* tmp;
	int d, n, i, j, ki = 0;
	points = kd_node->pfh;
	n = kd_node->n;
	d = 125;
	/* partition key index is that along which descriptors have most variance */
	for (j = 0; j < d; ++j)
	{
		mean = var = 0;
		for (i = 0; i < n; ++i)
			mean += points[i].fHistogram[j];
		mean +=
			mean /= n;
		for (i = 0; i < n; ++i)
		{
			x = points[i].fHistogram[j] - mean;
			var += x * x;
		}
		var /= n;
		if (var > var_max)
		{
			ki = j;
			var_max = var;
		}
	}
	/* partition key value is median of descriptor values at ki */
	tmp = (double*)calloc(n, sizeof(double));
	for (i = 0; i < n; ++i)
		tmp[i] = points[i].fHistogram[ki];
	kv = median_select(tmp, n);
	free(tmp);
	kd_node->ki = ki;
	kd_node->kv = kv;
}
/*
Finds the median value of an array.  The array's elements are re-ordered
by this function.

@param array an array; the order of its elememts is reordered
@param n number of elements in array

@return Returns the median value of array.
*/
double median_select(double* array, int n)
{
	return rank_select(array, n, (n - 1) / 2);
}
/*
Finds the element of a specified rank in an array using the linear time
median-of-medians algorithm by Blum, Floyd, Pratt, Rivest, and Tarjan.
The elements of the array are re-ordered by this function.

@param array an array; the order of its elelemts is reordered
@param n number of elements in array
@param r the zero-based rank of the element to be selected

@return Returns the element from array with zero-based rank r.
*/
double rank_select(double* array, int n, int r)
{
	double* tmp, med;
	int gr_5, gr_tot, rem_elts, i, j;

	/* base case */
	if (n == 1)
		return array[0];

	/* divide array into groups of 5 and sort them */
	gr_5 = n / 5;
	gr_tot = ceil(n / 5.0);
	rem_elts = n % 5;
	tmp = array;
	for (i = 0; i < gr_5; ++i)
	{
		insertion_sort(tmp, 5);
		tmp += 5;
	}
	insertion_sort(tmp, rem_elts);

	/* recursively find the median of the medians of the groups of 5 */
	tmp = (double*)calloc(gr_tot, sizeof(double));
	for (i = 0, j = 2; i < gr_5; ++i, j += 5)
		tmp[i] = array[j];
	if (rem_elts)
		tmp[i++] = array[n - 1 - rem_elts / 2];
	med = rank_select(tmp, i, (i - 1) / 2);
	free(tmp);

	/* partition around median of medians and recursively select if necessary */
	j = partition_array(array, n, med);
	if (r == j)
		return med;
	else if (r < j)
		return rank_select(array, j, r);
	else
	{
		array += j + 1;
		return rank_select(array, (n - j - 1), (r - j - 1));
	}
}
/*
Sorts an array in place into increasing order using insertion sort.

@param array an array
@param n number of elements
*/
void insertion_sort(double* array, int n)
{
	double k;
	int i, j;

	for (i = 1; i < n; ++i)
	{
		k = array[i];
		j = i - 1;
		while (j >= 0 && array[j] > k)
		{
			array[j + 1] = array[j];
			j -= 1;
		}
		array[j + 1] = k;
	}
}

/*
Partitions an array around a specified value.

@param array an array
@param n number of elements
@param pivot value around which to partition

@return Returns index of the pivot after partitioning
*/
int partition_array(double* array, int n, double pivot)
{
	double tmp;
	int p, i, j;

	i = -1;
	for (j = 0; j < n; ++j)
		if (array[j] <= pivot)
		{
			tmp = array[++i];
			array[i] = array[j];
			array[j] = tmp;
			if (array[i] == pivot)
				p = i;
		}
	array[p] = array[i];
	array[i] = pivot;

	return i;
}
void partition_nodes(struct KDNodePoint* kd_node)
{
	PFH* points;
	PFH tmp;
	double kv;
	int n, ki, p, i, j = -1;

	points = kd_node->pfh;
	n = kd_node->n;
	ki = kd_node->ki;
	kv = kd_node->kv;
	for (i = 0; i < n; ++i)
	{
		if (points[i].fHistogram[ki] <= kv)
		{
			tmp = points[++j];
			points[j] = points[i];
			points[i] = tmp;
			if (points[j].fHistogram[ki] == kv)
				p = j;
		}
	}
	tmp = points[p];
	points[p] = points[j];
	points[j] = tmp;

	/* if all records fall on same side of partition, make node a leaf */
	if (j == n - 1)
	{
		kd_node->leaf = 1;
		return;
	}

	kd_node->kd_left = kd_node_init(points, j + 1);
	kd_node->kd_right = kd_node_init(points + (j + 1), (n - j - 1));
}

void KDTreeRelease(KDNodePoint* kd_root)
{
	if (!kd_root)
		return;
	KDTreeRelease(kd_root->kd_left);
	KDTreeRelease(kd_root->kd_right);
	free(kd_root);
}
/*Finds an image feature's approximate k nearest neighbors in a kd tree using
Best Bin First search.

@param kd_root root of an image feature kd tree
@param point image feature for whose neighbors to search
@param k number of neighbors to find
@param nbrs pointer to an array in which to store pointers to neighbors
in order of increasing descriptor distance
@param max_nn_chks search is cut off after examining this many tree entries

@return Returns the number of neighbors found and stored in nbrs, or
-1 on error.*/
int KDTreeSearch(KDNodePoint* kd_root, PFH* point, int k, PFH*** nbrs, int max_nn_chks)
{
	KDNodePoint* expl;
	min_pq* min_pq;
	PFH* tree_point, ** _nbrs;
	bbf_data* bbf_data;
	int i, t = 0, n = 0;

	if (!nbrs || !point || !kd_root)
	{
		return -1;
	}

	_nbrs = (PFH**)calloc(k, sizeof(PFH*));
	min_pq = minpq_init();
	minpq_insert(min_pq, kd_root, 0);
	while (min_pq->n > 0 && t < max_nn_chks)
	{
		expl = (KDNodePoint*)minpq_extract_min(min_pq);
		if (!expl)
			goto fail;

		expl = explore_to_leaf(expl, point, min_pq);
		if (!expl)
			goto fail;
		for (i = 0; i < expl->n; ++i)
		{
			tree_point = &expl->pfh[i];
			bbf_data = (struct bbf_data*)malloc(sizeof(struct bbf_data));
			if (!bbf_data)
				goto fail;
			bbf_data->old_data = tree_point->feature_data;
			bbf_data->d = descr_dist_sq(point, tree_point);
			tree_point->feature_data = bbf_data;
			n += insert_into_nbr_array(tree_point, _nbrs, n, k);
		}
		t++;
	}

	minpq_release(&min_pq);
	for (i = 0; i < n; ++i)
	{
		bbf_data = (struct bbf_data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
		free(bbf_data);
	}
	*nbrs = _nbrs;
	return n;

fail:
	minpq_release(&min_pq);
	for (i = 0; i < n; ++i)
	{
		bbf_data = (struct bbf_data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
		free(bbf_data);
	}
	free(_nbrs);
	*nbrs = NULL;
	return -1;

}
/*
Explores a kd tree from a given node to a leaf.  Branching decisions are
made at each node based on the descriptor of a given feature.  Each node
examined but not explored is put into a priority queue to be explored
later, keyed based on the distance from its partition key value to the
given feature's desctiptor.

@param kd_node root of the subtree to be explored
@param point feature upon which branching decisions are based
@param min_pq a minimizing priority queue into which tree nodes are placed
as described above

@return Returns a pointer to the leaf node at which exploration ends or
NULL on error.
*/
struct KDNodePoint* explore_to_leaf(struct KDNodePoint* kd_node, struct PFH* point, struct min_pq* min_pq)
{
	struct KDNodePoint* unexpl;
	struct KDNodePoint* expl = kd_node;
	double kv;
	int ki;

	while (expl && !expl->leaf)
	{
		ki = expl->ki;
		kv = expl->kv;

		if (ki >= 125)
			return NULL;

		if (point->fHistogram[ki] <= kv)
		{
			unexpl = expl->kd_right;
			expl = expl->kd_left;
		}
		else
		{
			unexpl = expl->kd_left;
			expl = expl->kd_right;
		}

		if (minpq_insert(min_pq, unexpl, (int)abs(kv - point->fHistogram[ki])))
			return NULL;
	}
	return expl;
}
//直方图距离平方和
double descr_dist_sq(struct PFH* f1, struct PFH* f2)
{
	float diff, dsq = 0;
	float* descr1, * descr2;
	int i, d;

	d = 125;

	descr1 = f1->fHistogram;
	descr2 = f2->fHistogram;

	for (i = 0; i < d; ++i)
	{
		diff = descr1[i] - descr2[i];
		dsq += diff * diff;
	}
	return dsq;
}
/*
Inserts a feature into the nearest-neighbor array so that the array remains
in order of increasing descriptor distance from the search feature.

@param points feature to be inderted into the array; it's feature_data field
should be a pointer to a bbf_data with d equal to the squared descriptor
distance between feat and the search feature
@param nbrs array of nearest neighbors neighbors
@param n number of elements already in nbrs and
@param k maximum number of elements in nbrs

@return If feat was successfully inserted into nbrs, returns 1; otherwise
returns 0.
*/
int insert_into_nbr_array(struct PFH* point, struct PFH** nbrs, int n, int k)
{
	struct bbf_data* fdata, * ndata;
	double dn, df;
	int i, ret = 0;

	if (n == 0)
	{
		nbrs[0] = point;
		return 1;
	}

	/* check at end of array */
	fdata = (struct bbf_data*)point->feature_data;
	df = fdata->d;
	ndata = (struct bbf_data*)nbrs[n - 1]->feature_data;
	dn = ndata->d;
	if (df >= dn)
	{
		if (n == k)
		{
			point->feature_data = fdata->old_data;
			free(fdata);
			return 0;
		}
		nbrs[n] = point;
		return 1;
	}

	/* find the right place in the array */
	if (n < k)
	{
		nbrs[n] = nbrs[n - 1];
		ret = 1;
	}
	else
	{
		nbrs[n - 1]->feature_data = ndata->old_data;
		free(ndata);
	}
	i = n - 2;
	while (i >= 0)
	{
		ndata = (struct bbf_data*)nbrs[i]->feature_data;
		dn = ndata->d;
		if (dn <= df)
			break;
		nbrs[i + 1] = nbrs[i];
		i--;
	}
	i++;
	nbrs[i] = point;

	return ret;
}