#include "mac_grid.h"
#include "open_gl_headers.h" 
#include "camera.h"
#include "custom_output.h" 
#include "constants.h" 
#include <math.h>
#include <map>
#include <stdio.h>
#include <cstdlib>
#undef max
#undef min 
#include <fstream> 


// Globals
MACGrid target;

#define CHECK_DIVERGENCE


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = RenderMode::SHEETS;
bool MACGrid::theDisplayVel = false;//true

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 





MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);

   calculateAMatrix();
   calculatePreconditioner(AMatrix);
}

void MACGrid::initialize()
{
   reset();
}

double MACGrid::getDeltaTime() {
    double umax = 0;
	FOR_EACH_CELL {
                const double umag = getVelocity(getCenter(i,j,k)).Length();
                umax = umag > umax ? umag : umax;
	}
    //can use last iters max values since there's little delta between frames
    const double CFLfactor = 1.0;
    umax += sqrt(CFLfactor*theCellSize*(fbuoymax+fconfmax));
	return (CFLfactor*theCellSize / umax);
}

void MACGrid::updateSources()
{
    // Set initial values for density, temperature, velocity


    ////////////////first block
	int depth = theDim[MACGrid::Z] >> 1;
	double tempScale = 1.0;
    int width = 5;
    int startx = 6, stopx = startx + width, starty = 0, stopy = starty+width;
    for(int i=startx; i<stopx;i++){
        for(int j=starty; j<stopy; j++){
            mV(i,j+1,depth) = 6.0;
            mD(i,j,depth) = 1.0;
            mT(i,j,depth) = theBuoyancyAmbientTemperature*tempScale;
        }
    }

	//// Refresh particles in source.
	//for(int i=startx; i<stopx; i++) {
	//	for (int j = starty; j < stopy; j++) {
	//		for (int k = 0; k <= 0; k++) {
	//			vec3 cell_center(theCellSize*(i+0.5), theCellSize*(j+0.5), theCellSize*(k+0.5));
	//			for(int p=0; p<10; p++) {
    //                double a = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
    //                double b = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
    //                double c = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
    //                vec3 shift(a, b, c);
    //                vec3 xp = cell_center + shift;
    //                rendering_particles.push_back(xp);
    //            }
	//		}
	//	}
	//}

    //////////////////second block
    depth = theDim[MACGrid::Z] >> 1;
    tempScale = 1.4;
    startx = 27, stopx = startx + width, starty = 0, stopy = starty+width;
    for(int i=startx; i<stopx;i++){
        for(int j=starty; j<stopy; j++){
            mV(i,j+1,depth) = 4.0;
            mD(i,j,depth) = 1.0;
            mT(i,j,depth) = theBuoyancyAmbientTemperature*tempScale;
        }
    }

    //// Refresh particles in source.
    //for(int i=startx; i<stopx; i++) {
    //    for (int j = starty; j < stopy; j++) {
    //        for (int k = 0; k <= 0; k++) {
    //            vec3 cell_center(theCellSize*(i+0.5), theCellSize*(j+0.5), theCellSize*(k+0.5));
    //            for(int p=0; p<10; p++) {
    //                double a = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
    //                double b = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
    //                double c = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
    //                vec3 shift(a, b, c);
    //                vec3 xp = cell_center + shift;
    //                rendering_particles.push_back(xp);
    //            }
    //        }
    //    }
    //}
}

////////////////// BEGIN STEP 1/////////////////////

void MACGrid::advectVelocity(double dt)
{
    //Calculate new velocities and store in target

    //Get rid of these three lines after you implement yours
	//target.mU = mU;
    //target.mV = mV;
    //target.mW = mW;

    //Your code is here. It builds target.mU, target.mV and target.mW for all faces
	FOR_EACH_FACE {
				if (isValidFace(MACGrid::X, i, j, k)) {
					const double u = getVelocityX(getRewoundPosition(getFacePosition(MACGrid::X, i, j, k), dt));
					target.mU(i, j, k) = u;
				} if (isValidFace(MACGrid::Y, i, j, k)) {
					const double v = getVelocityY(getRewoundPosition(getFacePosition(MACGrid::Y, i, j, k), dt));
					target.mV(i, j, k) = v;
				} if (isValidFace(MACGrid::Z, i, j, k)) {
					const double w = getVelocityZ(getRewoundPosition(getFacePosition(MACGrid::Z, i, j, k), dt));
					target.mW(i, j, k) = w;
				}
	}


    // Then save the result to our object
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::addExternalForces(double dt) {
	computeBouyancy(dt);
	computeVorticityConfinement(dt);
}

void MACGrid::computeBouyancy(double dt)
{
	//Calculate bouyancy and store in target

	//Get rid of this line after you implement yours
	//target.mV = mV;

	//Your code is here. It modifies target.mV for all y face velocities.
	// f_buoy = vec3(0, -alpha*s + beta*(T - T_amb), 0);
	// Where alpha and beta are non-negative constants which can be set by
	// the user to achieve different behavior. Note that f_buoy = 0 where there's no
	// smoke and the temp is at ambient level.

	// copy the current state of mV as we need to keep accumulating
	target.mV = mV;

    //calc the force at grid centers then the face force is the ave of the 2 neighboring cells
	GridData fbuoy_c = GridData(); fbuoy_c.initialize(0.0);
	FOR_EACH_CELL {
				fbuoy_c(i,j,k) = -theBuoyancyAlpha * mD(i,j,k) +
								  theBuoyancyBeta * (mT(i,j,k) - theBuoyancyAmbientTemperature);
	}

	fbuoymax = 0;
    FOR_EACH_FACE {
				if(isValidFace(MACGrid::Y, i, j, k) && j != 0  && j != theDim[MACGrid::Y]) {
					//density and temp are defined at grid centers, will need to grab the two corresponding
					//cell centers forces, then average them and mult by dt (see fedkiw appendix page7 topleft), although..
					//is the interp done by getDensity and getTemp good enough

					//div around 90, seems more stable, weird artifacts at border
					//const vec3 facePos = getFacePosition(MACGrid::Y, i, j, k);
					//const double f_buoy = -theBuoyancyAlpha * getDensity(facePos) +
					//					  theBuoyancyBeta * (getTemperature(facePos) - theBuoyancyAmbientTemperature);
					//target.mV(i, j, k) += f_buoy; //dt not needed? looks better without at least...

					//using cell center forces, div around 90, seems less stable, but no weird artifacts at border
					const double fbuoy = 0.5 * (fbuoy_c(i,j,k) + fbuoy_c(i,j-1,k));
					fbuoymax = fbuoy > fbuoymax ? fbuoy : fbuoymax;
					target.mV(i, j, k) += dt*fbuoy;
				}
	}
	// and then save the result to our object
	mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	//Calculate vorticity confinement forces

	// Apply the forces to the current velocity and store the result in target
	// STARTED.

	//Get rid of this line after you implement yours
	//target.mU = mU;
	//target.mV = mV;
	//target.mW = mW;

	// Your code is here. It modifies target.mU,mV,mW for all faces.
    //copy the current state of vel components as we need to keep accumulating
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	//see fedkiw page 3 left, eq 9,10,11 and surrounding paras for vorticity confinement
	//see fedkiw appendix page 6 bottom right for vorticity discretization calc
	//see bridson siggraph page 58 and 59 for equations and discussion on vorticity confinement
    //w = curl of velocity field ("vorticity")
    //n = gradient of the magnitude of this vorticity
	//N = normalize(n)
	//f_conf = e * h * (N x w), e > 0 and controls confinement, h is cell size
	//calculate the force on grid centers then for each face take average of the 2 corresponding grid cell centers
	//border face can extrapolate or take the nearest grid center

	//need cell centered velocities for getting vorticity
	GridData mUc = GridData(); mUc.initialize(0.0);
	GridData mVc = GridData(); mVc.initialize(0.0);
	GridData mWc = GridData(); mWc.initialize(0.0);
	FOR_EACH_CELL {
                //faces are left biased(-1/2 corresponds to i,j,k and +1/2 corresponds to (x dir) i+1,j,k
				mUc(i,j,k) = 0.5 * (target.mU(i, j, k) + target.mU(i+1, j  , k  ));
				mVc(i,j,k) = 0.5 * (target.mV(i, j, k) + target.mV(i  , j+1, k  ));
				mWc(i,j,k) = 0.5 * (target.mW(i, j, k) + target.mW(i  , j  , k+1));
	}

    //fedkiw notation in appendixA for vorticity
	GridData w1 = GridData(); w1.initialize(0.0);
	GridData w2 = GridData(); w2.initialize(0.0);
	GridData w3 = GridData(); w3.initialize(0.0);
	GridData wMag = GridData(); wMag.initialize(0.0);
	const double invTwoCellSize = 1.0 / (2.0 * theCellSize);
	const double invCellSize = 1.0 / (theCellSize);
	FOR_EACH_CELL {
				const vec3 vort(
                        //the way fedkiw mentions, but produces divergent field
						//(mWc(i,j+1,k) - mWc(i,j-1,k) - mVc(i,j,k+1) + mVc(i,j,k-1)) * invTwoCellSize,
						//(mUc(i,j,k+1) - mUc(i,j,k-1) - mWc(i+1,j,k) + mWc(i-1,j,k)) * invTwoCellSize,
						//(mVc(i+1,j,k) - mVc(i-1,j,k) - mUc(i,j+1,k) + mUc(i,j-1,k)) * invTwoCellSize
						//same indexing except using cell faces, biased left, most divergence stability
						(target.mW(i,j+1,k) - target.mW(i,j-1,k) - target.mV(i,j,k+1) + target.mV(i,j,k-1)) * invTwoCellSize,
						(target.mU(i,j,k+1) - target.mU(i,j,k-1) - target.mW(i+1,j,k) + target.mW(i-1,j,k)) * invTwoCellSize,
						(target.mV(i+1,j,k) - target.mV(i-1,j,k) - target.mU(i,j+1,k) + target.mU(i,j-1,k)) * invTwoCellSize
                        //indexing for faces on either side of cell, no bias, less stable
						//(target.mW(i,j+1,k) - target.mW(i,j,k) - target.mV(i,j,k+1) + target.mV(i,j,k)) * invCellSize,
						//(target.mU(i,j,k+1) - target.mU(i,j,k) - target.mW(i+1,j,k) + target.mW(i,j,k)) * invCellSize,
						//(target.mV(i+1,j,k) - target.mV(i,j,k) - target.mU(i,j+1,k) + target.mU(i,j,k)) * invCellSize
				);
				w1(i,j,k) = vort[0];
				w2(i,j,k) = vort[1];
				w3(i,j,k) = vort[2];
				wMag(i,j,k) = vort.Length();
	}

    fconfmax = 0;
	FOR_EACH_CELL {
                //prefetch for later
				const vec3 vort(w1(i,j,k), w2(i,j,k), w3(i,j,k));

				const vec3 vortGrad(
						(wMag(i+1,j,k) - wMag(i-1,j,k)) * invTwoCellSize,
						(wMag(i,j+1,k) - wMag(i,j-1,k)) * invTwoCellSize,
						(wMag(i,j,k+1) - wMag(i,j,k-1)) * invTwoCellSize
				);
				const vec3 N = vortGrad / (vortGrad.Length() + 0.0000001);
                const vec3 fconf = theVorticityEpsilon * theCellSize * N.Cross(vort);

				//just push info to faces of the cell here instead of doing another loop
				if(isValidFace(MACGrid::X, i,j,k)   && i != 0 && i != theDim[MACGrid::X]) { target.mU(i,j,k)   += (dt*0.5*fconf[0]);}
				if(isValidFace(MACGrid::X, i+1,j,k) && i != 0 && i != theDim[MACGrid::X]) { target.mU(i+1,j,k) += (dt*0.5*fconf[0]);}
				if(isValidFace(MACGrid::Y, i,j,k)   && j != 0 && j != theDim[MACGrid::Y]) { target.mV(i,j,k)   += (dt*0.5*fconf[1]);}
				if(isValidFace(MACGrid::Y, i,j+1,k) && j != 0 && j != theDim[MACGrid::Y]) { target.mV(i,j+1,k) += (dt*0.5*fconf[1]);}
				if(isValidFace(MACGrid::Z, i,j,k)   && k != 0 && k != theDim[MACGrid::Z]) { target.mW(i,j,k)   += (dt*0.5*fconf[2]);}
				if(isValidFace(MACGrid::Z, i,j,k+1) && k != 0 && k != theDim[MACGrid::Z]) { target.mW(i,j,k+1) += (dt*0.5*fconf[2]);}
	}

	// Then save the result to our object
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;


//#ifdef CHECK_DIVERGENCE
//	GridData d = GridData(); d.initialize(0.0);
//	const double AMatrixCoeffDivide = (theAirDensity*theCellSize*theCellSize) / dt;
//	const double AMatrixCoeffDivideAndNegInvCellSize = AMatrixCoeffDivide * (-1.0 / theCellSize);
//	FOR_EACH_CELL {
//				d(i,j,k) = AMatrixCoeffDivideAndNegInvCellSize * ( (mU(i+1,j,k) - mU(i,j,k)) +
//																   (mV(i,j+1,k) - mV(i,j,k)) +
//																   (mW(i,j,k+1) - mW(i,j,k))
//				);
//			}
//
//	double total = 0;
//	FOR_EACH_CELL { total += d(i,j,k); }
//	if(total > 0.0000001) { PRINT_LINE("Non-zero Divergence: " << total); }
//#endif

}

void MACGrid::project(double dt)
{
	// Solve Ap = d for pressure
	// 1. Contruct d
	// 2. Construct A
	// 3. Solve for p
	// Subtract pressure from our velocity and save in target
	// STARTED.

	// Get rid of these 3 lines after you implement yours
	//target.mU = mU;
	//target.mV = mV;
	//target.mW = mW;

	// Your code is here. It solves for a pressure field and modifies target.mU,mV,mW for all faces.


	//construct d, see eq 4.22 on page 42 of bridsons course notes, see para 3 and 4 for explanation of AMatrix.
	// and the coefficients moved to the rhs of eq 4.22
	// velocities on the solid wall boundaries should be 0 before we set d
	// copy state of vel to target
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	GridData d = GridData(); d.initialize(0.0);
	const double AMatrixCoeffDivide = (theAirDensity*theCellSize*theCellSize) / dt;
	const double AMatrixCoeffDivideAndNegInvCellSize = AMatrixCoeffDivide * (-1.0 / theCellSize);
    FOR_EACH_CELL {
				d(i,j,k) = AMatrixCoeffDivideAndNegInvCellSize * ( (target.mU(i+1,j,k) - target.mU(i,j,k)) +
																   (target.mV(i,j+1,k) - target.mV(i,j,k)) +
											                       (target.mW(i,j,k+1) - target.mW(i,j,k))
											                     );
	}


	//construct A
	// A is constructed for solid static bounding box already
	// if there are any solid moving objects, AMatrix needs to be constructed
	// and its preconditioner calculated every time
	//calculateAMatrix();
	//calculatePreconditioner(AMatrix);


	//solve for p
    double tol = 0.00001;
	while(!preconditionedConjugateGradient(AMatrix, target.mP, d, 100, tol)) {tol *= 2; }
    //PRINT_LINE( "PCG converged with tolerance: " << tol );

	//subtract off the pressure gradients from the velocity fields.
	//bottom page 39 in bridson course notes, note: the face velocities for
	//solid boundaries must equal the velocity of teh solid boundary, in this case, 0
	const double c = dt / (theAirDensity * theCellSize);
	FOR_EACH_FACE {
				// set to vel of solid if face touches solid and solve for pressure there on other side
				// as described by 4.10 on page 39 of bridson course notes
				if(isValidFace(MACGrid::X, i,j,k)) {
					if (i == 0 || i == theDim[MACGrid::X]) {
						target.mU(i,j,k) = 0;
                    } else {
						target.mU(i,j,k) -= (c*(target.mP(i,j,k) - target.mP(i-1,j,k)));
					}
                }
				if(isValidFace(MACGrid::Y, i,j,k)) {
					if (j == 0 || j == theDim[MACGrid::Y]) {
						target.mV(i,j,k) = 0;
					} else {
						target.mV(i,j,k) -= (c*(target.mP(i,j,k) - target.mP(i,j-1,k)));
					}
				}
				if(isValidFace(MACGrid::Z, i,j,k)) {
					if (k == 0 || k == theDim[MACGrid::Z]) {
						target.mW(i,j,k) = 0;
					} else {
						target.mW(i,j,k) -= (c*(target.mP(i,j,k) - target.mP(i,j,k-1)));
					}
				}
	}

#ifdef CHECK_DIVERGENCE
	// Check border velocities:
	FOR_EACH_FACE {
				if (isValidFace(MACGrid::X, i, j, k)) {

					if (i == 0) {
						if (abs(target.mU(i,j,k)) > 0.0000001) {
							PRINT_LINE( "LOW X:  " << target.mU(i,j,k) );
							//target.mU(i,j,k) = 0;
						}
					}

					if (i == theDim[MACGrid::X]) {
						if (abs(target.mU(i,j,k)) > 0.0000001) {
							PRINT_LINE( "HIGH X: " << target.mU(i,j,k) );
							//target.mU(i,j,k) = 0;
						}
					}

				}
				if (isValidFace(MACGrid::Y, i, j, k)) {


					if (j == 0) {
						if (abs(target.mV(i,j,k)) > 0.0000001) {
							PRINT_LINE( "LOW Y:  " << target.mV(i,j,k) );
							//target.mV(i,j,k) = 0;
						}
					}

					if (j == theDim[MACGrid::Y]) {
						if (abs(target.mV(i,j,k)) > 0.0000001) {
							PRINT_LINE( "HIGH Y: " << target.mV(i,j,k) );
							//target.mV(i,j,k) = 0;
						}
					}

				}
				if (isValidFace(MACGrid::Z, i, j, k)) {

					if (k == 0) {
						if (abs(target.mW(i,j,k)) > 0.0000001) {
							PRINT_LINE( "LOW Z:  " << target.mW(i,j,k) );
							//target.mW(i,j,k) = 0;
						}
					}

					if (k == theDim[MACGrid::Z]) {
						if (abs(target.mW(i,j,k)) > 0.0000001) {
							PRINT_LINE( "HIGH Z: " << target.mW(i,j,k) );
							//target.mW(i,j,k) = 0;
						}
					}
				}
			}
	// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
	// TODO: Fix duplicate code:
	FOR_EACH_CELL {
				// Construct the vector of divergences d:
				double velLowX = target.mU(i,j,k);
				double velHighX = target.mU(i+1,j,k);
				double velLowY = target.mV(i,j,k);
				double velHighY = target.mV(i,j+1,k);
				double velLowZ = target.mW(i,j,k);
				double velHighZ = target.mW(i,j,k+1);
				double divergence = ((velHighX - velLowX) + (velHighY - velLowY) + (velHighZ - velLowZ)) / theCellSize;
				if (abs(divergence) > 0.02 ) {
					PRINT_LINE("WARNING: Divergent! ");
					PRINT_LINE("Divergence: " << divergence);
					PRINT_LINE("Cell: " << i << ", " << j << ", " << k);
				}
			}
#endif


	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;


}
////////////////// END STEP 1/////////////////////

void MACGrid::advectTemperature(double dt)
{
    //Calculate new temp and store in target

    //Get rid of this line after you implement yours
    //target.mT = mT;

    //Your code is here. It builds target.mT for all cells.
    FOR_EACH_CELL {
		target.mT(i,j,k) = getTemperature( getRewoundPosition( getCenter(i,j,k), dt) );
	}

    // Then save the result to our object
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    //Calculate new densitities and store in target

    //Get rid of this line after you implement yours
    //target.mD = mD;

    //Your code is here. It builds target.mD for all cells.
	FOR_EACH_CELL {
		target.mD(i,j,k) = getDensity( getRewoundPosition( getCenter(i,j,k), dt) );
	}

    // Then save the result to our object
    mD = target.mD;
}


void MACGrid::advectRenderingParticles(double dt) {
	rendering_particles_vel.resize(rendering_particles.size());
	for (size_t p = 0; p < rendering_particles.size(); p++) {
		vec3 currentPosition = rendering_particles[p];
        vec3 currentVelocity = getVelocity(currentPosition);
        vec3 nextPosition = currentPosition + currentVelocity * dt;
        vec3 clippedNextPosition = clipToGrid(nextPosition, currentPosition);
        // Keep going...
        vec3 nextVelocity = getVelocity(clippedNextPosition);
        vec3 averageVelocity = (currentVelocity + nextVelocity) / 2.0;
        vec3 betterNextPosition = currentPosition + averageVelocity * dt;
        vec3 clippedBetterNextPosition = clipToGrid(betterNextPosition, currentPosition);
        rendering_particles[p] = clippedBetterNextPosition;
		rendering_particles_vel[p] = averageVelocity;
	}
}


vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}


vec3 MACGrid::getRewoundPosition(const vec3 & currentPosition, const double dt) {

	/*
	// EULER (RK1):
	vec3 currentVelocity = getVelocity(currentPosition);
	vec3 rewoundPosition = currentPosition - currentVelocity * dt;
	vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
	return clippedRewoundPosition;
	*/

	// HEUN / MODIFIED EULER (RK2):
	vec3 currentVelocity = getVelocity(currentPosition);
	vec3 rewoundPosition = currentPosition - currentVelocity * dt;
	vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
	// Keep going...
	vec3 rewoundVelocity = getVelocity(clippedRewoundPosition);
	vec3 averageVelocity = (currentVelocity + rewoundVelocity) / 2.0;
	vec3 betterRewoundPosition = currentPosition - averageVelocity * dt;
	vec3 clippedBetterRewoundPosition = clipToGrid(betterRewoundPosition, currentPosition);
	return clippedBetterRewoundPosition;

}


vec3 MACGrid::clipToGrid(const vec3& outsidePoint, const vec3& insidePoint) {
	/*
	// OLD:
	vec3 rewindPosition = outsidePoint;
	if (rewindPosition[0] < 0) rewindPosition[0] = 0; // TEMP!
	if (rewindPosition[1] < 0) rewindPosition[1] = 0; // TEMP!
	if (rewindPosition[2] < 0) rewindPosition[2] = 0; // TEMP!
	if (rewindPosition[0] > theDim[MACGrid::X]) rewindPosition[0] = theDim[MACGrid::X]; // TEMP!
	if (rewindPosition[1] > theDim[MACGrid::Y]) rewindPosition[1] = theDim[MACGrid::Y]; // TEMP!
	if (rewindPosition[2] > theDim[MACGrid::Z]) rewindPosition[2] = theDim[MACGrid::Z]; // TEMP!
	return rewindPosition;
	*/

	vec3 clippedPoint = outsidePoint;

	for (int i = 0; i < 3; i++) {
		if (clippedPoint[i] < 0) {
			vec3 distance = clippedPoint - insidePoint;
			double newDistanceI = 0 - insidePoint[i];
			double ratio = newDistanceI / distance[i];
			clippedPoint = insidePoint + distance * ratio;
		}
		if (clippedPoint[i] > getSize(i)) {
			vec3 distance = clippedPoint - insidePoint;
			double newDistanceI = getSize(i) - insidePoint[i];
			double ratio = newDistanceI / distance[i];
			clippedPoint = insidePoint + distance * ratio;
		}
	}

#ifdef _DEBUG
	// Make sure the point is now in the grid:
	if (clippedPoint[0] < 0 || clippedPoint[1] < 0 || clippedPoint[2] < 0 || clippedPoint[0] > getSize(0) || clippedPoint[1] > getSize(1) || clippedPoint[2] > getSize(2)) {
		PRINT_LINE("WARNING: Clipped point is outside grid!");
	}
#endif

	return clippedPoint;

}


double MACGrid::getSize(int dimension) {
	return theDim[dimension] * theCellSize;
}


int MACGrid::getCellIndex(int i, int j, int k)
{
	return i + j * theDim[MACGrid::X] + k * theDim[MACGrid::Y] * theDim[MACGrid::X];
}


int MACGrid::getNumberOfCells()
{
	return theDim[MACGrid::X] * theDim[MACGrid::Y] * theDim[MACGrid::Z];
}


bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}


bool MACGrid::isValidFace(int dimension, int i, int j, int k)
{
	if (dimension == 0) {
		if (i > theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
			return false;
		}
	} else if (dimension == 1) {
		if (i >= theDim[MACGrid::X] || j > theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
			return false;
		}
	} else if (dimension == 2) {
		if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k > theDim[MACGrid::Z]) {
			return false;
		}
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}



vec3 MACGrid::getFacePosition(int dimension, int i, int j, int k)
{
	if (dimension == 0) {
		return vec3(i * theCellSize, (j + 0.5) * theCellSize, (k + 0.5) * theCellSize);
	} else if (dimension == 1) {
		return vec3((i + 0.5) * theCellSize, j * theCellSize, (k + 0.5) * theCellSize);
	} else if (dimension == 2) {
		return vec3((i + 0.5) * theCellSize, (j + 0.5) * theCellSize, k * theCellSize);
	}

	return vec3(0,0,0); //???

}

void MACGrid::calculateAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

bool MACGrid::preconditionedConjugateGradient( const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	// Initial guess p = 0
	FOR_EACH_CELL { p(i, j, k) = 0.0; }

	GridData r = d; // Residual vector.

	/*
	PRINT_LINE("r: ");
	FOR_EACH_CELL {
		PRINT_LINE(r(i,j,k));
	}
	*/
	GridData z; z.initialize();
	applyPreconditioner(r, A, z); // Auxillary vector.
	/*
	PRINT_LINE("z: ");
	FOR_EACH_CELL {
		PRINT_LINE(z(i,j,k));
	}
	*/

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma; // According to TA. Here???

		apply(A, s, z); // z = applyA(s);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);
		//p += alpha * s;

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);
		//r -= alpha * z;

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations, tol: " << tolerance );
			return true; //return p;
		}

		applyPreconditioner(r, A, z); // z = applyPreconditioner(r);

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	//PRINT_LINE( "PCG didn't converge with tolerance: " << tolerance );
	return false;
}


void MACGrid::calculatePreconditioner(const GridDataMatrix & A) {

	precon.initialize();

    // Build the modified incomplete Cholesky preconditioner following Fig 4.2 on page 36 of Bridson's 2007 SIGGRAPH fluid course notes.
    //       This corresponds to filling in precon(i,j,k) for all cells
	const double tau = 0.97;
	FOR_EACH_CELL {
				double one = A.plusI(i-1, j, k) * precon(i-1,j,k);
				double two = A.plusJ(i, j-1, k) * precon(i,j-1,k);
				double thr = A.plusK(i, j, k-1) * precon(i,j,k-1);
                one *= one;
				two *= two;
				thr *= thr;

				const double four = A.plusI(i-1,j,k)*(A.plusJ(i-1,j,k)+A.plusK(i-1,j,k))*precon(i-1,j,k)*precon(i-1,j,k);
				const double five = A.plusJ(i,j-1,k)*(A.plusI(i,j-1,k)+A.plusK(i,j-1,k))*precon(i,j-1,k)*precon(i,j-1,k);
				const double six  = A.plusK(i,j,k-1)*(A.plusI(i,j,k-1)+A.plusJ(i,j,k-1))*precon(i,j,k-1)*precon(i,j,k-1);

				const double e = A.diag(i,j,k) - one - two - thr - tau*(four+five+six);

				precon(i,j,k) = 1.0 / sqrt(e + 0.0000001);
	}
}


void MACGrid::applyPreconditioner(const GridData & r, const GridDataMatrix & A, GridData & z) {

    // change if(0) to if(1) after you implement calculatePreconditoner function.

    if(1) {

        // APPLY THE PRECONDITIONER:
        // Solve Lq = r for q:
        GridData q;
        q.initialize();
        FOR_EACH_CELL {
                    //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
                    double t = r(i, j, k) - A.plusI(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k)
                               - A.plusJ(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k)
                               - A.plusK(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                    q(i, j, k) = t * precon(i, j, k);
                    //}
                }
        // Solve L^Tz = q for z:
        FOR_EACH_CELL_REVERSE {
                    //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
                    double t = q(i, j, k) - A.plusI(i, j, k) * precon(i, j, k) * z(i + 1, j, k)
                               - A.plusJ(i, j, k) * precon(i, j, k) * z(i, j + 1, k)
                               - A.plusK(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                    z(i, j, k) = t * precon(i, j, k);
                    //}
                }
    }
    else{
        // Unpreconditioned CG: Bypass preconditioner:
        z = r;
        return;
    }

}



double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}


void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}


void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}


void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}


double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}


void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}

void MACGrid::saveParticle(std::string filename){
	Partio::ParticlesDataMutable *parts = Partio::create();
	Partio::ParticleAttribute posH, vH;
	posH = parts->addAttribute("position", Partio::VECTOR, 3);
	vH = parts->addAttribute("v", Partio::VECTOR, 3);
	for (unsigned int i = 0; i < rendering_particles.size(); i++)
	{
		int idx = parts->addParticle();
		float *p = parts->dataWrite<float>(posH, idx);
		float *v = parts->dataWrite<float>(vH, idx);
		for (int k = 0; k < 3; k++)
		{
			p[k] = rendering_particles[i][k];
			v[k] = rendering_particles_vel[i][k];
		}
	}
	
	Partio::write(filename.c_str(), *parts);
	parts->release();
}

void MACGrid::saveDensity(std::string filename){
	Partio::ParticlesDataMutable *density_field = Partio::create();
	Partio::ParticleAttribute posH, rhoH;
	posH = density_field->addAttribute("position", Partio::VECTOR, 3);
	rhoH = density_field->addAttribute("density", Partio::VECTOR, 1);
	FOR_EACH_CELL{
		int idx = density_field->addParticle();
		float *p = density_field->dataWrite<float>(posH, idx);
		float *rho = density_field->dataWrite<float>(rhoH, idx);
		vec3 cellCenter = getCenter(i, j, k);
		for (int l = 0; l < 3; l++)
		{
			p[l] = cellCenter[l];
		}
		rho[0] = getDensity(cellCenter);
	}
	Partio::write(filename.c_str(), *density_field);
	density_field->release();
}

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0);
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
           //vel.Normalize(); 
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{
	
	double value = mD(i, j, k); 
	vec4 coldColor(0.5, 0.5, 1.0, value);
	vec4 hotColor(1.0, 0.5, 0.5, value);
    return LERP(coldColor, hotColor, mT(i, j, k));
	

	/*
	// OLD:
    double value = mD(i, j, k); 
    return vec4(1.0, 0.9, 1.0, value);
	*/
}

vec4 MACGrid::getRenderColor(const vec3& pt)
{
	double value = getDensity(pt);
	vec4 coldColor(0.5, 0.5, 1.0, value);
	vec4 hotColor(1.0, 0.5, 0.5, value);
    return LERP(coldColor, hotColor, getTemperature(pt));

	/*
	// OLD:
    double value = getDensity(pt); 
    return vec4(1.0, 1.0, 1.0, value);
	*/
}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}


void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}