// PairMatchAdapter.cpp : 定义控制台应用程序的入口点。
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <pcl/io/pcd_io.h>
#include <iostream>
#include <string.h>
#include "KDTree.h"
#include <pcl/common/transforms.h>
#include <pcl/registration/transformation_estimation_svd.h>

using namespace pcl;
using namespace std;

//TargetCloud:根据SourceCloud用KDTree搜索排序后的点云特征点；输出刚体旋转矩阵
registration::TransformationEstimationSVD<PointXYZ, PointXYZ, float>::Matrix4 KDTreeRansac(PointCloud<PointXYZ>::Ptr TargetCloudFeature, PointCloud<PointXYZ>::Ptr SourceCloudFeature)
{
	//确定nCorrs组匹配点对
	const int nCorrs = 3;
	PointCloud<PointXYZ>::Ptr TargetPointsChosen(new PointCloud<PointXYZ>);
	PointCloud<PointXYZ>::Ptr SourcePointsChosen(new PointCloud<PointXYZ>);
	TargetPointsChosen->width = nCorrs;
	TargetPointsChosen->height = 1;
	TargetPointsChosen->is_dense = true;
	TargetPointsChosen->resize(TargetPointsChosen->width * TargetPointsChosen->height);
	SourcePointsChosen->width = nCorrs;
	SourcePointsChosen->height = 1;
	SourcePointsChosen->is_dense = true;
	SourcePointsChosen->resize(SourcePointsChosen->width * SourcePointsChosen->height);
	int nPointsIdx[nCorrs];//点对索引
	int nRand = 0;
	double MinFitnessScore = 1000000.0;
	registration::TransformationEstimationSVD<PointXYZ, PointXYZ>::Matrix4 FinalTransformMatrix;//最终旋转矩阵
	PointCloud<PointXYZ>::Ptr SourceCloudFeatureTransformed(new PointCloud<PointXYZ>());
	srand((unsigned int)time(0));
	for (int nIter = 0; nIter < 2000; nIter++)
	{
		ShowProgress(nIter, 2000);
		//cout << "iter:" << nIter << "------------" << endl;
		for (int i = 0; i < nCorrs; i++)
		{
			nRand = rand();
			int p = nRand % (TargetCloudFeature->size());
			nPointsIdx[i] = p;
		}
		for (int i = 0; i < nCorrs; ++i)
		{
			TargetPointsChosen->points[i] = TargetCloudFeature->at(nPointsIdx[i]);
			SourcePointsChosen->points[i] = SourceCloudFeature->at(nPointsIdx[i]);
		}
		//根据点对生成旋转矩阵
		registration::TransformationEstimationSVD<PointXYZ, PointXYZ> TransEstimSVD;
		registration::TransformationEstimationSVD<PointXYZ, PointXYZ, float>::Matrix4 TransformMatrix;
		TransEstimSVD.estimateRigidTransformation(*SourcePointsChosen, *TargetPointsChosen, TransformMatrix);
		/*cout << "estimated rotation:" << endl;
		printf("    | %6.3f %6.3f %6.3f | \n", TransformMatrix(0, 0), TransformMatrix(0, 1), TransformMatrix(0, 2));
		printf("R = | %6.3f %6.3f %6.3f | \n", TransformMatrix(1, 0), TransformMatrix(1, 1), TransformMatrix(1, 2));
		printf("    | %6.3f %6.3f %6.3f | \n", TransformMatrix(2, 0), TransformMatrix(2, 1), TransformMatrix(2, 2));
		printf("t = < %0.3f, %0.3f, %0.3f >\n", TransformMatrix(0, 3), TransformMatrix(1, 3), TransformMatrix(2, 3));*/
		//计算匹配分数
		double dFitnessScore = 0.0;
		transformPointCloud(*SourceCloudFeature, *SourceCloudFeatureTransformed, TransformMatrix);
		KdTreeFLANN<PointXYZ> kdtree;
		kdtree.setInputCloud(TargetCloudFeature);
		int K = 1;//找到最邻近的K个点
		std::vector<int> vNearestNeiborIdxs(K);
		std::vector<float> vNearestNeiborDists(K);
		int nNearestNeiborNum = 0;
		for (int i = 0; i < SourceCloudFeatureTransformed->points.size(); i = i + 1)
		{
			kdtree.nearestKSearch(SourceCloudFeatureTransformed->points[i], 1, vNearestNeiborIdxs, vNearestNeiborDists);
			if (vNearestNeiborDists[0] <= 100000000)
			{
				dFitnessScore += vNearestNeiborDists[0];
				nNearestNeiborNum++;
			}
		}
		if (nNearestNeiborNum > 0)
			dFitnessScore = (dFitnessScore / nNearestNeiborNum);
		else
			dFitnessScore = std::numeric_limits<double>::max();
		//cout << "FitnessScore is:" << dFitnessScore << endl;
		if (dFitnessScore < MinFitnessScore)
		{
			MinFitnessScore = dFitnessScore;
			FinalTransformMatrix = TransformMatrix;
		}
	}
	return FinalTransformMatrix;//返回最终旋转矩阵
}

int main(int argc, char* argv[])
{
	//argv[1] = "e:\\OfficeEF\\Middle\\OfficeE77.poi";
	//argv[2] = "e:\\OfficeEF\\Middle\\OfficeF20.poi";
	//argv[3] = "e:\\OfficeEF\\Middle\\OfficeE77_OfficeF20.tie";
	if (argv[1] == NULL || argv[2] == NULL || argv[3] == NULL)
	{
		cout << "insufficient parameters.";
		return 0;
	}
	//读取文件内poi信息到点云的PFH结构体数组里并查询大小
	int* nTargetNum = new int, * nSourceNum = new int;
	int* nTargetID = new int;
	int* nSourceID = new int;
	PFH* TargetPFH = LoadPFHFile(argv[1], nTargetNum, nTargetID);
	PFH* SourcePFH = LoadPFHFile(argv[2], nSourceNum, nSourceID);
	if (TargetPFH == NULL || SourcePFH == NULL)
	{
		cout << "failed to open poi file.";
		return 0;
	}
	//建立KD搜索树
	KDNodePoint* RootTargetPFH;
	RootTargetPFH = KDTreeBuild(TargetPFH, *nTargetNum);
	PFH** NeiborPFHs;
	PointCloud<PointXYZ>::Ptr TargetFeature(new PointCloud<PointXYZ>);
	PointCloud<PointXYZ>::Ptr SourceFeature(new PointCloud<PointXYZ>);
	PointXYZ point;
	double dDist;//点对PFH距离
	cout << "calculating " << argv[1] << " with " << argv[2] << "...";
	for (int i = 0; i < *nSourceNum; ++i)
	{
		KDTreeSearch(RootTargetPFH, &SourcePFH[i], 1, &NeiborPFHs, 200);
		dDist = 0.0;
		for (int j = 0; j < 125; j++)
			dDist += (SourcePFH[i].fHistogram[j] - NeiborPFHs[0]->fHistogram[j]) * (SourcePFH[i].fHistogram[j] - NeiborPFHs[0]->fHistogram[j]);
		if (dDist < 20.0 && dDist != 0)//判断特征点特征空间距离是否足够小
		{
			point.x = NeiborPFHs[0]->x; point.y = NeiborPFHs[0]->y; point.z = NeiborPFHs[0]->z;
			TargetFeature->push_back(point);
			point.x = SourcePFH[i].x; point.y = SourcePFH[i].y; point.z = SourcePFH[i].z;
			SourceFeature->push_back(point);
		}
	}
	registration::TransformationEstimationSVD<PointXYZ, PointXYZ, float>::Matrix4 FinalTransMatrix = KDTreeRansac(TargetFeature, SourceFeature);//最终旋转矩阵	

	PointCloud<PointXYZ>::Ptr SourceFeatureTransed(new PointCloud<PointXYZ>);
	transformPointCloud(*SourceFeature, *SourceFeatureTransed, FinalTransMatrix);
	//visualize_pcd(TargetFeature, SourceFeatureTransed, SourceFeatureTransed, SourceFeatureTransed);
	//根据最终旋转矩阵剔除错误对应点
	PointCloud<PointXYZ>::Ptr TargetCorrs(new PointCloud<PointXYZ>);//正确的对应点
	PointCloud<PointXYZ>::Ptr SourceCorrs(new PointCloud<PointXYZ>);
	double dPointDis = 0.0;
	for (int i = 0; i < TargetFeature->size(); ++i)
	{
		dPointDis = (TargetFeature->at(i).x - SourceFeatureTransed->at(i).x) * (TargetFeature->at(i).x - SourceFeatureTransed->at(i).x)
			+ (TargetFeature->at(i).y - SourceFeatureTransed->at(i).y) * (TargetFeature->at(i).y - SourceFeatureTransed->at(i).y)
			+ (TargetFeature->at(i).z - SourceFeatureTransed->at(i).z) * (TargetFeature->at(i).z - SourceFeatureTransed->at(i).z);
		if (dPointDis < 400)
		{
			TargetCorrs->push_back(TargetFeature->at(i));
			SourceCorrs->push_back(SourceFeature->at(i));
		}
	}
	//保存对应点文件tie
	FILE* pfTie;
	fopen_s(&pfTie, argv[3], "w");
	if (pfTie == NULL)
	{
		cout << "failed to open tie file.";
		return 0;
	}
	fprintf(pfTie, "%d %d %d\n", TargetCorrs->size(), *nTargetID, *nSourceID);
	for (int i = 0; i < TargetCorrs->size(); ++i)
		fprintf(pfTie, "%.4f %.4f %.4f   %.4f %.4f %.4f\n",
			TargetCorrs->at(i).x / 1000.0, TargetCorrs->at(i).y / 1000.0, TargetCorrs->at(i).z / 1000.0, SourceCorrs->at(i).x / 1000.0, SourceCorrs->at(i).y / 1000.0, SourceCorrs->at(i).z / 1000.0);
	//清理 
	fclose(pfTie);
	cout << "success." << endl;
	delete[] nTargetNum;
	nTargetNum = NULL;
	delete[] nSourceNum;
	nSourceNum = NULL;
	delete[] nTargetID;
	nTargetID = NULL;
	delete[] nSourceID;
	nSourceID = NULL;
	KDTreeRelease(RootTargetPFH);
	return 0;
}