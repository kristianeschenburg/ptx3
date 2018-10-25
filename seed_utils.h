/*
 * seed_utils.h
 *
 *  Created on: Oct 22, 2018
 *      Author: keschenb
 */

#ifndef SEED_UTILS_H_
#define SEED_UTILS_H_

#include <iostream>
#include <math.h>
#include <fstream>
#include "newimage/newimageall.h"
#include "meshclass/meshclass.h"
#include "probtrackxOptions.h"
#include "spherical.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace std;

namespace seed_utilities {

int round(float coordinate) {

	int rounded = (int) MISCMATHS::round(coordinate);
	return rounded;
}

int voxel_value(volume<short int> cortex, ColumnVector point) {
	return cortex(round(point(1)), round(point(2)), round(point(3)));
}

// Convert spherical coordinates to Cartesian coordinates
ColumnVector sphere2cart(float theta, float phi) {

	ColumnVector cart(3);
	float x,y,z;

	x = sin(theta)*cos(phi);
	y = sin(theta)*sin(phi);
	z = cos(theta);

	cart << x << y << z;

	return cart;

}

// Convert Cartesian coordinates to spherical coordinates
ColumnVector cart2sphere(ColumnVector cartesian) {

	float x, y, z, theta, phi;
	ColumnVector spherical(2);
	x = cartesian(1);
	y = cartesian(2);
	z = cartesian(3);

	theta = acos(z);
	phi = atan2(y,x);

	spherical << theta << phi;

	return spherical;

}

ColumnVector scale_angle(ColumnVector angle, float steplength, volume<short int> cortex) {

	ColumnVector scaled(3);
	float x,y,z;
	x = angle(1);
	y = angle(2);
	z = angle(3);

	float xscale, yscale, zscale;
	xscale = x*steplength / cortex.xdim();
	yscale = y*steplength / cortex.ydim();
	zscale = z*steplength / cortex.zdim();

	scaled << xscale << yscale << zscale;

	return scaled;

}

float euclidean(ColumnVector source, ColumnVector target) {

	float euc;
	ColumnVector diff = target - source;

	euc = sqrt(pow(diff(1),2.0) + pow(diff(2),2.0) + pow(diff(3),2.0));

	return euc;
}

volume<float> step_over_volume(ColumnVector seed, ColumnVector angle,
		volume<float> streamline_counts) {

	ColumnVector upd(3);

	upd << seed(1) << seed(2) << seed(3);
	streamline_counts(round(upd(1)),round(upd(2)),round(upd(3))) = 1;
	upd = upd + 2*angle;
	streamline_counts(round(upd(1)),round(upd(2)),round(upd(3))) = 3;
	upd = upd + 2*angle;
	streamline_counts(round(upd(1)),round(upd(2)),round(upd(3))) = 5;

	return streamline_counts;
}

}

#endif /* SEED_UTILS_H_ */
