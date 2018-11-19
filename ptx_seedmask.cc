/*  Copyright (C) 2004 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk

    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK


    LICENCE

    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")

    The Software remains the property of the University of Oxford ("the
    University").

    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.

    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.

    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.

    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    Innovation@innovation.ox.ac.uk quoting reference DE/9564. */

#include "miscmaths/miscmaths.h"
#include "streamlines.h"
#include "ptx_seedmask.h"
#include "seed_utils.h"
#include <time.h>

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;
using namespace MISCMATHS;
using namespace seed_utilities;


void seedmask()
{ 

	probtrackxOptions& opts =probtrackxOptions::getInstance();

	// we need a reference volume for CSV
	// (in case seeds are a list of surfaces)
	volume<short int> refvol;
	if(opts.seedref.value()!="")
		read_volume(refvol,opts.seedref.value());
	else
		read_volume(refvol,opts.maskfile.value());

	if (opts.save_subpaths.value()) {
		cout << "Save SubPaths Value: " << opts.save_subpaths.value() << endl;
	}



	cout<<"load seeds"<<endl;
	CSV seeds(refvol);
	seeds.set_convention(opts.meshspace.value());
	seeds.load_rois(opts.seedfile.value());

	cout<<"done\n"<<endl;
	if(seeds.nVols()==0 && opts.seedref.value()==""){
		cerr<<"Warning: need to set a reference volume when defining a surface-based seed"<<endl;
	}


	// KE
	// If initial directions file supplied,
	// load angles and initialize inner and outer surfaces
	Matrix directions;
	SeedUtilities SU;
	CSV outer(refvol), inner(refvol);
	if (opts.forceangle.value()) {

		cout << "load initial direction" << endl;
		directions = read_ascii_matrix(opts.initdir.value());
		cout << "# Directions: " << directions.Nrows() << endl;
		cout << "# angle columns: " << directions.Ncols() << endl;
		cout << "done\n" << endl;

		cout << "load inner and outer surface" << endl;

		outer.set_convention(opts.meshspace.value());
		outer.load_rois(opts.outersurf.value());

		inner.set_convention(opts.meshspace.value());
		inner.load_rois(opts.innersurf.value());

		cout << "done\n" << endl;
	}



	cout << "FlipSign Value: " << opts.flipsign.value() << endl;


	// Initialize tracking classes
	Streamliner  stline      (seeds);
	Counter      counter     (stline);
	counter.initialise();
	Seedmanager  seedmanager (counter);

	/* counter.initialise() initializes the output matrices
	 * that store the streamline counts
	*/

	srand(opts.rseed.value()); // need to reinitialise random seed because of GIFTI!!


	// Initialize an empty volume of same size as seedref
	// This volume will save the initial streamline direction in the first step
	volume<float> direction_counts;
	direction_counts.reinitialize(stline.get_seeds().xsize(),
			stline.get_seeds().ysize(), stline.get_seeds().zsize());
	copybasicproperties(stline.get_seeds().get_refvol(), direction_counts);
	direction_counts = 0;
	//
	//
	//

	int keeptotal=0;

	time_t _time;
	_time=time(NULL);
	// seed from volume-like ROIs
	if(seeds.nVols()>0){
		cout << "Volume seeds" << endl;
		vector<int> triangles; //to avoid connections between vertices of the same triangle...but not used for volumes
		triangles.push_back(-1);

		for(int roi=1;roi<=seeds.nVols();roi++){
			cout<<"volume "<<roi<<endl;

			for(int z=0;z<seeds.zsize();z++){
				if(opts.verbose.value()>=1)
					cout <<"sl "<<z<<endl;
				for(int y=0;y<seeds.ysize();y++){
					for(int x=0;x<seeds.xsize();x++){
						if(seeds.isInRoi(x,y,z,roi)){
							counter.updateSeedLocation(seeds.get_volloc(roi-1,x,y,z),-1,triangles);
							if(opts.verbose.value()>=1){
								cout <<"run"<<endl;
								cout <<x<<" "<<y<<" "<<z<<endl;
							}
							keeptotal += seedmanager.run((float)x,(float)y,(float)z,
									false,-1,opts.sampvox.value());
						}
					}
				}
			}
		}
	}

	cout << "OneWayOnly flag: " << opts.onewayonly.value() << endl;

	// seed from surface-like ROIs
	if(seeds.nSurfs()>0){
		cout << "Surface seeds" << endl;
		ColumnVector pos;
		ColumnVector pos_gw;
		for(int i=0;i<seeds.nSurfs();i++){
			cout<<"surface "<<i<<endl;

			// inform user if whole surface is used or not
			if( seeds.nActVertices(i) != seeds.nVertices(i) ){
				cout << "  Using a subset of the vertices labeled active (i.e. non zero value)" << endl;
				cout << "   set all values to 0 or non-zero to use entire surface" << endl;
			}

			// KE
			// Main loop here over each seed point
			// For each seed p, generate opts.nparticles.value() streamlines

			for(int p=0;p<seeds.get_mesh(i).nvertices();p++){

				// check if active point
				if(seeds.get_mesh(i).get_pvalue(p)==0.0)
					continue;

				//to avoid connections between vertices of the same triangle
				CsvMpoint vertex=seeds.get_mesh(i).get_point(p);
				vector<int> triangles;
				for(int t=0;t<vertex.ntriangles();t++){
					triangles.push_back(vertex.get_trID(t));
				}

				counter.updateSeedLocation(seeds.get_surfloc(i,p), i, triangles);
				pos = seeds.get_vertex_as_vox(i,p);

				// KE
				// if forcing tracking with initial directions
				if (opts.forceangle.value()) {

					// KE
					// get spherical coordinates, convert to Cartesian
					float theta = directions(p+1,1);
					float phi = directions(p+1,2);

					ColumnVector angle = SU.sphere2cart(theta, phi);
					angle = SU.normalize(angle);

					// get inner and outer surface positions
					ColumnVector inner_pos = inner.get_vertex_as_vox(0, p);
					ColumnVector outer_pos = outer.get_vertex_as_vox(0, p);

					// KE
					// Distance from G/W to Pial (cortical thickness)
					// Distance from seed to G/W
					// Distance from seed to Pial
					float thickness = SU.euclidean(inner_pos,outer_pos);
					float to_outers = SU.euclidean(pos, outer_pos);

					// KE
					if (to_outers < thickness) {
						pos = inner_pos;
						to_outers = thickness;
					}

					// KE
					ColumnVector move_pos = pos + angle;

					float upd_outers = SU.euclidean(move_pos, outer_pos);
					if (upd_outers < to_outers) {
						angle = (-1)*angle;
					}

					// KE
					// convert back to spherical coordinates
					angle = SU.normalize(angle);
					ColumnVector spherical = SU.cart2sphere(angle);
					theta = spherical(1);
					phi = spherical(2);

					// KE
					// set streamliner object with initial angle
					seedmanager.get_stline().set_angles(theta,phi);
					/*
					 * End KE
					 */

				}

				if(opts.verbose.value()>=1){
					cout <<"run"<<endl;
					cout << pos(1) << " " << pos(2) << " " << pos(3) << endl;
				}

				// perform tracking
				keeptotal += seedmanager.run(pos(1), pos(2), pos(3),
						opts.onewayonly.value(), -1, opts.sampvox.value());

			}
		}
	}

	cout << endl << "time spent tracking: " << (time(NULL)-_time) << " seconds"<< endl << endl;

	// save results
	cout << "save results" << endl;
	counter.save_total(keeptotal);
	counter.save();

	cout<<"finished"<<endl;
}


