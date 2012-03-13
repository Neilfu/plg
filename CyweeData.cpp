// CyweeData.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#define N 50000
#define T 100  // smooth
#define M 5   // 5轴手柄
#define GravityThreshold 10
#define FEATURE_NUM_1 8
#define THRESHOLD 0.1


int flag[N],cnt = 0,border[N][2],bcnt;
double dis[N];

bool smooth(int k,int f)
{
	int i,m;
	double n1 = 0,n2 = 0;
	if (k+T > cnt)
		m = cnt;
	else
		m = k+T;

	for (i = k+1;i < m;i++)
	{
		if (flag[i] == f)
			n1 ++;
		else
			n2 ++;
	}
	if(n1/n2 >= 2)
		return true;
	return false;
}

void segment()
{
	ifstream fin;
	fin.open("RawData.txt");
	ofstream fout;
	fout.open("test.txt");

	double G = 269,x,y,z;
	string str;
	int i = 0;
	bool begin  = false;

	getline(fin,str);
	while (fin >> x)
	{
		x -= 2048;
		fin >> y;
		y -= 2048;
		fin >> z;
		z -= 2048;
		dis[cnt] = sqrt(x*x + y*y + z*z);
		//fout << x << ' ' << y << ' ' << z << ' ' << dis[cnt] << endl;
		getline(fin,str);
		if (abs(G-dis[cnt]) < GravityThreshold)
			flag[cnt] = 0;
		else
			flag[cnt] = 1;
		cnt ++;
	}
	while (i < cnt)
	{
		if (!begin)
		{
			if (flag[i] == 1)
			{
				if (smooth(i,1))
				{
					begin = true;
					border[bcnt][0] = i;
				}
			}
		} 
		else
		{
			if (flag[i] == 0)
			{
				if (smooth(i,0))
				{
					begin = false;
					border[bcnt][1] = i;
					bcnt ++;
				}
			}
		}
		i ++;
	}
	cout << bcnt << endl;
	fin.close();
	fout.close();
}

void ExtractFeature_1()
{
	ofstream fout;
	fout.open("feature1.txt");
	double ** feature;	
	int i,j,k;
	feature = new double * [bcnt];
	for (i = 0;i < bcnt;i ++)
	{
		feature[i] = new double[FEATURE_NUM_1];
	}
	for (i = 0;i < bcnt;i++)
	{
		k = (border[i][1] - border[i][0] + 1) / FEATURE_NUM_1;
		for (j = 0;j < FEATURE_NUM_1;j ++)
		{
			feature[i][j] = dis[border[i][0] + j * k];
			fout << feature[i][j] << ' ';
		}
		fout << endl;
	}
	fout.close();
}

void ExtractParam(int m,double *t,double &compress_ratio,double &avg_inc,double &avg_dec,double &avg_dis,double &var)
{
	int i = border[m][0],k,section_len,section_num = 0,start = border[m][0],cdis = 0;
	double ax,ay,Lxx,Lxy,Q,F,a,b,a0,data_len;
	avg_inc = 0;
	avg_dec = 0;
	avg_dis = 0;
	while (i <= border[m][1])
	{
		section_num ++;
		section_len = 1;
		Lxx = 0;
		Lxy = 0;
		Q = 0;
		a0 = 0;
		ax = t[i];
		ay = dis[i];
		avg_dis += dis[i];cdis ++;
		F = sqrt(Q/section_len) / ax;
		i ++;
		while (i <= border[m][1])
		{
			ax = section_len * ax / (section_len + 1) + t[i] / (section_len + 1);
			ay = section_len * ay / (section_len + 1) + dis[i] / (section_len + 1);
			Lxx += section_len * (ax - t[i]) * (ax - t[i]) / (section_len + 1);
			Lxy += section_len * (ax - t[i]) * (ay - dis[i]) / (section_len + 1);
			section_len ++;			
			a = Lxy / Lxx; 
			b = ay - a * ax;
			Q = 0;
			for (k = start;k < start + section_len;k ++)
			{
				Q += pow((a * t[k] + b - dis[k]),2);
			}
			F = sqrt(Q/section_len) / ax;
			if (F > THRESHOLD)
			{
				start = i;
				if (a0 > 0)
					avg_inc += a0;
				else
					avg_dec += a0;
				break;
			}
			avg_dis += dis[i]; cdis ++;
			a0 = a;
			i ++;
		}
	}
	data_len = border[m][1] - border[m][0] + 1;
	compress_ratio = data_len / section_num;
	avg_inc = avg_inc / compress_ratio;
	avg_dec = avg_dec / compress_ratio;
	avg_dis = avg_dis / data_len;
	var = 0;
	for (i = border[m][0];i <= border[m][1];i ++)
	{
		var += pow(dis[i]-avg_dis,2);
	}
	var = var / pow(data_len,2);
}

void ExtractFeature_2()
{
	int i,j;
	double compress_ratio,avg_inc,avg_dec,avg_dis,var,t[N];
	ofstream fout;
	fout.open("feature2.txt");
	for (i = 0;i < bcnt;i ++)
	{
		for (j = border[i][0];j <= border[i][1];j ++)
		{
			t[j] = j - border[i][0] + 1;
		}
		ExtractParam(i,t,compress_ratio,avg_inc,avg_dec,avg_dis,var);
		fout << compress_ratio << ' ' << avg_inc << ' ' << avg_dec << ' ' << avg_dis << ' ' << var << endl;
	}
	fout.close();
}

int _tmain(int argc, _TCHAR* argv[])
{
	system("WaveTest.exe");
	segment();
	for (int i = 0;i < bcnt;i ++)
	{
		cout << i + 1 << ' ' << border[i][0] << ' ' << border[i][1] << ' ' << border[i][1] - border[i][0] << endl;
	}
	ExtractFeature_1();
	ExtractFeature_2();

	system("pause");
	return 0;
}
