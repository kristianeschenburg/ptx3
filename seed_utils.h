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
#include <cmath>


namespace seed_utilities {
class SeedUtilities {

	public:

		/*
		 * Round float to nearest int
		 */
		int round(float coordinate) {

			int rounded = (int) MISCMATHS::round(coordinate);
			return rounded;
		}

		/*
		 * Return value of volume at voxel
		 */
		int voxel_value(volume<short int> cortex, ColumnVector point) {
			return cortex(round(point(1)), round(point(2)), round(point(3)));
		}

		/*
		 * Compute vector length
		 */
		float vector_length(ColumnVector cart_pos) {

			float length = sqrt(pow(cart_pos(1), 2.0) +
								pow(cart_pos(2), 2.0) +
								pow(cart_pos(3), 2.0));
			return length;
		}

		/*
		 * Normalize vector to unit length
		 */
		ColumnVector normalize(ColumnVector cart_pos) {

			float length = vector_length(cart_pos);
			ColumnVector normed = cart_pos/length;

			return normed;

		}

		// Convert spherical coordinates to Cartesian coordinates
		ColumnVector sphere2cart(ColumnVector spheres) {

			ColumnVector cart(3);
			float x,y,z;
			float theta, phi;
			theta = spheres(1);
			phi = spheres(2);

			x = sin(theta)*cos(phi);
			y = sin(theta)*sin(phi);
			z = cos(theta);

			cart << x << y << z;
			cart = normalize(cart);

			return cart;

		}

		// Convert Cartesian coordinates to spherical coordinates
		ColumnVector cart2sphere(ColumnVector cart_pos) {

			float x, y, z, theta, phi;
			ColumnVector spherical(2);

			cart_pos = normalize(cart_pos);
			x = cart_pos(1);
			y = cart_pos(2);
			z = cart_pos(3);

			theta = acos(z);
			phi = atan2(y,x);

			spherical << theta << phi;

			return spherical;

		}

		/*
		 * Scale direction vector by volume voxel dimensions
		 */
		ColumnVector scale_angle(ColumnVector angle, float steplength, volume<short int> cortex) {

			ColumnVector scaled(3);
			float x,y,z;
			angle = normalize(angle);

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

		/*
		 * Compute Euclidean distance between two vectors
		 */
		float euclidean(ColumnVector source, ColumnVector target) {

			float euc;
			ColumnVector diff = target - source;
			euc = vector_length(diff);

			return euc;
		}

		/*
		 * Check if seed position is closer to pial surface than
		 * corresponding g/w matter seed.  If it is, seed pos
		 * becomes g/w seed.
		 *
		 * Parameters:
		 * 		seed_pos: current seed position
		 * 		inner_pos: analogous seed position on g/w surface
		 * 		outer_pos: analogous seed position on pial surface
		 */
		ColumnVector check_seed(ColumnVector seed_pos, ColumnVector inner_pos,
				ColumnVector outer_pos) {

			float thickness = euclidean(inner_pos, outer_pos);
			float to_outers = euclidean(seed_pos, outer_pos);

			if (to_outers < thickness) {
				seed_pos = inner_pos;
				to_outers = thickness;
			}

			return seed_pos;

		}

		/*
		 * Check that the angle is in the proper direction
		 * 	i.e. inwards towards the center of the brain
		 * 	Parameters:
		 * 		seed_pos: current seed position
		 * 		outer_pos: analogous seed position on pial surface
		 * 		angle: direction in which to move
		 */
		ColumnVector check_angle(ColumnVector seed_pos, ColumnVector outer_pos,
				ColumnVector angle) {

			float to_outers = euclidean(seed_pos, outer_pos);
			ColumnVector move_pos = seed_pos + angle;

			float upd_outers = euclidean(move_pos, outer_pos);
			if (upd_outers < to_outers) {
				angle = (-1)*angle;
			}
			return angle;
		}

		/*
		 * Given a direction vector, take steps in that direction
		 * Assign each voxel that this movement crosses a value
		 * Return this stepped-over volume
		 *
		 * Parameters:
		 * 		seed: current seed position
		 * 		angle: direction to move in
		 * 		streamline_counts: volume to update with steps
		 */
		volume<float> step_over_volume(ColumnVector seed, ColumnVector angle,
				volume<float> streamline_counts) {

			ColumnVector upd(3);

			upd << seed(1) << seed(2) << seed(3);
			for (int i = 0; i < 5; i++) {
				streamline_counts(round(upd(1)),round(upd(2)),round(upd(3))) = i+1;
				upd = upd + angle;
			}
			return streamline_counts;
		}

		/*
		 * Compute the great circle distance between two
		 * vectors of spherical coordinates.
		 *
		 * Assumes:
		 * 		sphere(1) = azimuthal angle
		 * 		sphere(2) = polar angle
		 *
		 * 		that the spherical angles were computed from
		 * 		normalized unit vectors
		 *
		 * Parameters:
		 * 		sph_1, sph_2: 2 spherical coordinates
		 */
		float spherical_distance(ColumnVector sph_1, ColumnVector sph_2) {

			/*
			 * 	cos_inv [
			 * 			cos(theta_1)*cos(theta_2) +
			 * 			sin(theta_1)*sin(theta_2)
			 * 			cos(phi_1-phi_2) ]
			 */

			float inner = cos(sph_1(1))*cos(sph_2(1)) +
					sin(sph_1(1))*sin(sph_2(1)) *
					cos(sph_1(2)-sph_2(2));

			float great_circle = acos(inner);

			return great_circle;
		}
};
}



#endif /* SEED_UTILS_H_ */
