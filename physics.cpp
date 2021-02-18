/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world* jello, struct point a[8][8][8])
{
	/* for you to implement ... */
	// F = ma

	// Hook's law: double f_hook=kx
	// rest length: structural 1.0/7, shear: 1.0 / 7 * sqrt(2) OR 1.0 / 7 * sqrt(3), bend: 1.0 / 7 * 2
	static const double restLengthStructural = 1.0 / 7;
	static const double restLengthShearDiagonal2D = restLengthStructural * sqrt(2.0);
	static const double restLengthShearDiagonal3D = restLengthStructural * sqrt(3.0);
	static const double restLengthBend = restLengthStructural * 2;

	// for both linear Hook and linear damp, the forces on the spring are always opposite on the two ends
	// so we need to update both: cache[x1][y1][z1][x2][y2][z2] == -cache[x2][y2][z2][x1][y1][z1] == the value
	static bool fCacheForceCalculated[8][8][8][8][8][8];
	static struct point fCacheForceBetweenPoints [8][8][8][8][8][8];

	// Linear Hook: double f_hookLinear = k_hook * x
	memset(fCacheForceCalculated, 0, sizeof(bool) * 512 * 512);
	memset(fCacheForceBetweenPoints, 0, sizeof(struct point) * 512 * 512);
	struct point fHookLinear[8][8][8];
	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				// hookLinear: structural
				// hookLinear: shear2D
				// hookLinear: shear3D
				// hookLinear: bend
			}
		}
	}

	// Linear damping: double f_linearDamp=-k_damp * v, where v is the relative speed ON THE SPRING'S DIRECTION
	memset(fCacheForceCalculated, 0, sizeof(bool) * 512 * 512);
	memset(fCacheForceBetweenPoints, 0, sizeof(struct point) * 512 * 512);
	struct point fDampLinear[8][8][8];
	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				// Damplinear: structural
				// Damplinear: shear2D
				// Damplinear: shear3D
				// Damplinear: bend
			}
		}
	}

	// Force field: double f_extForceField
	// trilinear interpolate and get the actual force at certain points
	//struct point fExtForce[8][8][8];

	// at last: double f_all = ma, then a = f_all / m
	
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world* jello)
{
	int i, j, k;
	point a[8][8][8];

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
				jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
				jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
				jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
				jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
				jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

			}
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world* jello)
{
	point F1p[8][8][8], F1v[8][8][8],
		F2p[8][8][8], F2v[8][8][8],
		F3p[8][8][8], F3v[8][8][8],
		F4p[8][8][8], F4v[8][8][8];

	point a[8][8][8];


	struct world buffer;

	int i, j, k;

	buffer = *jello; // make a copy of jello

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				pMULTIPLY(jello->v[i][j][k], jello->dt, F1p[i][j][k]);
				pMULTIPLY(a[i][j][k], jello->dt, F1v[i][j][k]);
				pMULTIPLY(F1p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F1v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F2p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F2p[i][j][k]);
				// F2v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F2v[i][j][k]);
				pMULTIPLY(F2p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F2v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F3p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F3v[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 1.0, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 1.0, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);


	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F4p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F4v[i][j][k]);

				pMULTIPLY(F2p[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1p[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4p[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->p[i][j][k], jello->p[i][j][k]);

				pMULTIPLY(F2v[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4v[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->v[i][j][k], jello->v[i][j][k]);
			}

	return;
}
