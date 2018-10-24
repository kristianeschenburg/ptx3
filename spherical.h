/*
 * spherical.h
 *
 *  Created on: Oct 16, 2018
 *      Author: keschenb
 */

#ifndef SPHERICAL_H_
#define SPHERICAL_H_

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

namespace cartsphere {

class Spherical {

	vector <vector<double> > angles;
	string coord_file;

	public:

	Spherical(string in_file): coord_file(in_file) {
		get_angles();
	}

	Spherical() { }
	~Spherical() { }

	void get_angles() {

		ifstream file (coord_file.c_str());
		string line;

		while (getline(file, line)) {
			string line_value;
			vector<double> line_values;
			stringstream ss(line);

			while (getline(ss, line_value, '\t')) {

				stringstream slv;
				double f = 0.0;

				slv << line_value;
				slv >> f;
				line_values.push_back(f);
			}
			angles.push_back(line_values);
		}
	}

	double theta(int sample) {
		return angles.at(sample).at(0);
	}
	double phi(int sample) {
		return angles.at(sample).at(1);
	}

	vector<double> point(int sample) {
		return angles.at(sample);
	}

	vector <vector <double> > samples() {
		return angles;
	}

};

}

#endif /* SPHERICAL_H_ */
