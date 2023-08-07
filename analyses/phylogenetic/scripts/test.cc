#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <regex>


template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	os << "[";
	for (auto const& x : v)
		os << x << " ";
	return os << "]";
}

#include "rangetree.h"


using namespace std;

vector<Point> range_points(int n, int m) {
	vector<Point> points;
	for (int i=0; i<n; i++)
		for (int j=0; j<m; j++)
			points.push_back(Point(i, j));
	return points;
}


void test1() {
	vector<Point> points = range_points(1, 10);
	Point p(1, 1);
	Segment s(points.begin(), points.end());
	assert(s.query(0, 10) == INF);
	s.modify(0, 2);
	assert(s.query(0, 10) == 2);
	s.modify(2, 6);
	assert(s.query(0, 10) == 2);
	s.modify(9, 1);
	assert(s.query(0, 10) == 1);
	cerr << "test 1 done" << endl;
}

void test1_1() {
	vector<Point> points {Point(1, 0), Point(1, 2), Point(1, 9)};
	Point p(1, 1);
	Segment s(points.begin(), points.end());
	assert(s.query(0, 10) == INF);
	s.modify(0, 2);
	assert(s.query(0, 10) == 2);
	s.modify(2, 6);
	assert(s.query(0, 10) == 2);
	s.modify(9, 1);
	assert(s.query(0, 10) == 1);
	cerr << "test 1_1 done" << endl;
}

int min_range(const vector<vector<int>>& a, int x1, int y1, int x2, int y2) {
	int r = INF;
	for (int i=x1; i<x2; i++)
		for (int j=y1; j<y2; j++)
			r = min(r, a[i][j]);
	return r;
}

void test2() {
	int n = 10, m = 10;
	vector<Point> points = range_points(n, m);
	cerr << "rt creating ... " << endl;
	RangeTree rt(n, m, points);
	cerr << "rt created " << endl;
	vector<vector<int>> naiv(n, vector<int>(m, INF));
	for (int t=0; t<1000; t++) {
		if (rand() % 10 == 0) {
			int x = rand() % n, y = rand() % m;
			naiv[x][y] = min(naiv[x][y], rand() % 100);
			rt.modify(x, y, naiv[x][y]);
		} else {
			int x1 = rand() % n, y1 = rand() % m;
			int x2 = rand() % n, y2 = rand() % m;
			if (min_range(naiv, x1, y1, x2, y2) != rt.query(x1, y1, x2, y2)) {
				cerr << "Q " << min_range(naiv, x1, y1, x2, y2) << " " << rt.query(x1, y1, x2, y2) << " [" << x1 << "," << y1 << "-" << x2 << "," << y2 << "]" << endl;
				for (int i=0; i<n; i++) {
					for  (int j=0; j<m; j++) {
						cerr << naiv[i][j] << " " ;
					}
					cerr << endl;
				}
			}
			assert(min_range(naiv, x1, y1, x2, y2) == rt.query(x1, y1, x2, y2));
		}
	}
}

void test2_2() {
	int n = 10, m = 10;
	vector<Point> points;
	for (int i=0; i<12; i++)
		points.push_back(Point( (i * 65537) % n , (i * 9601) % m));
	cerr << "rt creating ... " << points << endl;
	RangeTree rt(n, m, points);
	cerr << "rt created " << endl;
	vector<vector<int>> naiv(n, vector<int>(m, INF));
	for (int t=0; t<1000; t++) {
		if (rand() % 10 == 0) {
			//int x = rand() % n, y = rand() % m;
			int ind = rand() % points.size();
			int x = points[ind].x, y = points[ind].y;
			naiv[x][y] = min(naiv[x][y], rand() % 100);
			rt.modify(x, y, naiv[x][y]);
		} else {
			int x1 = rand() % n, y1 = rand() % m;
			int x2 = rand() % n, y2 = rand() % m;
			if (min_range(naiv, x1, y1, x2, y2) != rt.query(x1, y1, x2, y2)) {
				cerr << "Q " << min_range(naiv, x1, y1, x2, y2) << " " << rt.query(x1, y1, x2, y2) << " [" << x1 << "," << y1 << "-" << x2 << "," << y2 << "]" << endl;
				for (int i=0; i<n; i++) {
					for  (int j=0; j<m; j++) {
						cerr << naiv[i][j] << " " ;
					}
					cerr << endl;
				}
			}
			assert(min_range(naiv, x1, y1, x2, y2) == rt.query(x1, y1, x2, y2));
		}
	}
}

void test3() {
	int n = 10, m = 10;
	for (int i=0; i<n; i++)
		for (int j=0; j<m; j++) {
			vector<vector<int>> naiv(n, vector<int>(m, INF));
			vector<Point> points = range_points(n, m);
			RangeTree rt(n, m, points);
			rt.modify(i, j, 5);
			naiv[i][j] = 5;
			for (int x1=0; x1<n; x1++)
			for (int x2=x1; x2<n; x2++)
			for (int y1=0; y1<n; y1++)
			for (int y2=y1; y2<n; y2++) {
			assert(min_range(naiv, x1, y1, x2, y2) == rt.query(x1, y1, x2, y2));
			}
			//cout << rt.query(0, 0, 2, 2) << endl;
		}
}

void test4() {
	vector<vector<int>> naiv = {{2, 5, 17, 85}, 
		{INF, 31, 65, 86}, 
			{97, 36, INF, 64}, 
				{INF, 57, 54, INF}}; 

	int n = 4, m = 4;
	vector<Point> points = range_points(n, m);
	RangeTree rt(n, m, points);
	for (int i=0; i<n; i++) {
		for  (int j=0; j<m; j++) {
			rt.modify(i, j, naiv[i][j]);
		}
	}
	cout << rt.query(0, 1, 2, 2) << endl;

}

int days_since(const string& date, int y = 1900, int m = 1, int d = 1) {
        tm then = {0};
        then.tm_year = y - 1900;
        then.tm_mon = m - 1;
        then.tm_mday = d;
        time_t then_secs = mktime(&then);

        tm date_tm = {0};
        if (!strptime(date.c_str(), "%Y-%m-%d", &date_tm)) {
		throw exception();
	}
	cerr << "D " << date_tm.tm_year << " " << date_tm.tm_mon << " " << date_tm.tm_mday << " " << date_tm.tm_hour << " " << endl;
	cerr << "D " << then.tm_year << " " << then.tm_mon << " " << then.tm_mday << " " << then.tm_hour << " " << endl;
        time_t t = mktime(&date_tm);  // t is now your desired time_t
        
        double days = difftime(t, then_secs);
	cerr << "days_since_1900 " << date << " " << days << endl;
	return ((long)days) / 24 / 60 / 60;
}       

double years(const string& date) {
	tm date_tm = {0};
	if (!strptime(date.c_str(), "%Y-%m-%d", &date_tm)) {
		throw exception();
	}
	mktime(&date_tm);  // t is now your desired time_t

	return (double) date_tm.tm_yday / 366.0 + date_tm.tm_year + 1900;
}

string years_to_calendar(double years) {
	tm date_tm = {0};
	date_tm.tm_year = int(years) - 1900;
	date_tm.tm_mday = int( 366 * (years-int(years))) + 1;
	//cerr << date_tm.tm_mon << " " << date_tm.tm_mday << " " << date_tm.tm_yday << endl;
	mktime(&date_tm);  // t is now your desired time_t
	//cerr << date_tm.tm_mon << " " << date_tm.tm_mday << " " << date_tm.tm_yday << endl;
	char buff[100];
	strftime(buff, sizeof(buff), "%Y-%m-%d", &date_tm);
	return buff;
}

int month_since(const string& date, int y = 1900, int m = 1, int d = 1) {
	tm date_tm = {0};
	if (!strptime(date.c_str(), "%Y-%m-%d", &date_tm)) {
		throw exception();
	}
	return (date_tm.tm_year + 1900 - y) * 12 + (date_tm.tm_mon- m - 1);
}       

int main(void) {
	/*
	test1();
	test1_1();
	test2();
	test2_2();
	test3();
	test4();
	*/

	/*
	cout << days_since("1900-01-01") << endl;
	cout << days_since("2000-01-01") << endl;
	cout << days_since("1901-01-02") << endl;
	cout << days_since("2020-01-01") << endl;
	cout << days_since("2020-01-02") << endl;
	*/

	/*
	static const std::regex re{"^[0-9]{4}-[0-9]{2}-[0-9]{2}$"};
	cout << std::regex_match("2020-01-0", re) << endl;
	*/

	/*
	cout << "Y " << years("1900-01-01") << endl;
	cout << "Y " << years("2020-01-01") << endl;
	cout << "Y " << years("2020-06-01") << endl;
	cout << "Y " << years("2020-12-29") << endl;
	*/

	//cerr << years_to_calendar(2000) << endl;

	cerr << month_since("2020-01-01") << " " <<
		month_since("2020-01-30") << " " <<
		month_since("2020-02-01") << " " <<
		month_since("2021-01-01") << " " <<
		endl;

	return 0;
}
