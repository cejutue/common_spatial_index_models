#include "stdafx.h"
#include <vector>
#include <list>
#include <queue>
#include <map>
#include <algorithm>
#include <iostream>
#include <deque>
#include <math.h>
#include <time.h>
#include <functional>  
#include <memory>
#include <set>
#include <map>
#include <queue>
#include <string>
#include <memory>
#include <functional>
#include <math.h>
using namespace std;
typedef unsigned long long ull;

#include "orderingindex.cpp"




int testz2()
{
	Z2SpaceFillingCurve f(31);
	std::vector<double> coords = { 78,15,112,45 };
	ull g = f.index(78, 15);
	double x, y;
	//正反算测试
	f.invert(g, x, y);
	g = f.index(90, 19);
	f.invert(g, x, y);
	g = f.index(115, 54);
	f.invert(g, x, y);
	g = f.index(75, 19);
	f.invert(g, x, y);
	g = f.index(79, 13);
	f.invert(g, x, y);
	g = f.index(79, 46);
	f.invert(g, x, y);

	ZRange r1(f.index(78, 15), f.index(112, 45));
	ZRange rin(f.index(78, 15), f.index(112, 45));
	ZRange rover(f.index(90, 19), f.index(115, 54));
	ZRange rout(f.index(113, 15), f.index(118, 56));
	Z2 z2(2);
	//点测试
	bool tt = z2.contains(r1, f.index(90, 19));
	bool t1f = z2.contains(r1, f.index(115, 54));
	bool t2f = z2.contains(r1, f.index(75, 19));
	bool t3f = z2.contains(r1, f.index(79, 13));
	bool t4f = z2.contains(r1, f.index(79, 46));
	//矩形测试 
	bool true1 = z2.overlaps(r1, rin);
	bool true2 = z2.overlaps(r1, rover);
	bool false1 = z2.overlaps(r1, rout);

	//XZOrdering
	return 0;
}

int testXZ()
{
	XZ2SpaceFillingCurve f(31, -180, -90, 180, 90);
	ull g = f.index(78, 15, 112, 45);

	XZ2SpaceFillingCurve sfc(12, -180, -90, 180, 90);
	ull poly = sfc.index(10, 10, 12, 12);
	ull polyright = sfc.index(12.1, 10, 14, 12);
	ull polyleft = sfc.index(8, 10, 10, 12);
	ull polytop = sfc.index(10, 13, 12, 15);
	ull polybottom = sfc.index(10, 8, 12, 10);
	ull polybottomright = sfc.index(12, 8, 14, 10);

	std::vector<QueryWindow> containing_overlapping = {
		QueryWindow(9.0, 9.0, 13.0, 13.0),
		QueryWindow(-180.0, -90.0, 180.0, 90.0),
		QueryWindow(0.0, 0.0, 180.0, 90.0),
		QueryWindow(0.0, 0.0, 20.0, 20.0) ,
		QueryWindow(11.0, 11.0, 13.0, 13.0),
		QueryWindow(9.0, 9.0, 11.0, 11.0),
		QueryWindow(10.5, 10.5, 11.5, 11.5),
		QueryWindow(11.0, 11.0, 11.0, 11.0)
	};
	// note: in general, some disjoint ranges will match due to false positives
	std::vector<QueryWindow>  disjoint = {
		QueryWindow(-180.0, -90.0, 8.0, 8.0),
		QueryWindow(0.0, 0.0, 8.0, 8.0),
		QueryWindow(9.0, 9.0, 9.5, 9.5),
		QueryWindow(20.0, 20.0, 180.0, 90.0) };
	for (auto item  : containing_overlapping)
	{
		std::vector<QueryWindow> q;
		q.push_back(item);
		auto ranges = sfc.ranges(q);
		for (auto itemr  : ranges)
		{
			if (itemr.lower <= poly && itemr.upper <= poly)
			{
				printf("$bbox -  match");
			}
			else
			{
				printf("$bbox - no match");
			}
		}

	}
	//forall(disjoint) {
	//	bbox = >
	//		val ranges = sfc.ranges(Seq(bbox)).map(r = > (r.lower, r.upper))
	//		val matches = ranges.exists(r = > r._1 <= poly && r._2 >= poly)
	//		if (matches) {
	//			logger.warn(s"$bbox - invalid match")
	//		}
	//	matches must beFalse

	return 0;
}

int main()
{

	//testz2();
	//testXZ();
	return 0;
}