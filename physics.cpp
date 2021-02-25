/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

static bool isValidIndex(int i, int j, int k);
static void computeLinearHookAndDamp(struct world* jello, struct point fHookLinear[8][8][8], struct point fDampLinear[8][8][8]);

static struct point computeHookForceOnA(struct world* jello, int xA, int yA, int zA, int xB, int yB, int zB, double restLength,
	bool fCacheForceHookCalculated[8][8][8][8][8][8], struct point fCacheForceHookLinear[8][8][8][8][8][8]);
static struct point computeDampForceOnA(struct world* jello, int xA, int yA, int zA, int xB, int yB, int zB,
	bool fCacheForceDampCalculated[8][8][8][8][8][8], struct point fCacheForceDampLinear[8][8][8][8][8][8]);

void computeCollisionHookAndDampFromBoundingBox(struct world* jello, struct point fCollision[8][8][8]);
static struct point computeCollisionHookForce(const struct point& massPoint, double kCollision,
	double collisionPointX, double collisionPointY, double collisionPointZ);
// here we simplified the case, which is that the velocityB is always zero
// because for now all the objects other than the cube are static
static struct point computeCollisionDampForce(const struct point& massPoint, const struct point& curPointVelocity, double dCollision,
	double collisionPointX, double collisionPointY, double collisionPointZ);

void computeExternalForces(struct world* jello, struct point fExtForce[8][8][8]);
inline int getForceFieldIndex(int forceGridIndexI, int forceGridIndexJ, int forceGridIndexK, int resolution) {
	return (forceGridIndexI * resolution * resolution) + (forceGridIndexJ * resolution) + (forceGridIndexK);
}

void computeCollisionHookAndDampFromIncPlane(struct world* jello, struct point fCollision[8][8][8]);

/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world* jello, struct point a[8][8][8])
{
	// F = ma
	struct point fHookLinear[8][8][8] = { {{0.0}} };
	struct point fDampLinear[8][8][8] = { {{0.0}} };
	computeLinearHookAndDamp(jello, fHookLinear, fDampLinear);

	// collision: bounding box: hook & damp
	struct point fCollision[8][8][8] = { {{0.0}} };
	computeCollisionHookAndDampFromBoundingBox(jello, fCollision);

	// collision: inclined plane: F(x,y,z) = ax+by+cz+d = 0, on the plane. 
	// If more than 50% points is on one side and the rest is on the other side
	// then the minority of points are colliding with the plane
	// then we generate artificial collision springs for them
	if (jello->incPlanePresent) {
		computeCollisionHookAndDampFromIncPlane(jello, fCollision);
	}

	// Force field: double f_extForceField
	// trilinear interpolate and get the actual force at certain points
	struct point fExtForce[8][8][8] = { {{0.0}} };
	// the force field is between [-2, 2]^3
	computeExternalForces(jello, fExtForce);

	// at last: double f_all = ma, then a = f_all / m
	double dragForceMultiplierSinglePoint;
	if (jello->integrator[0] == 'E')
	{
		dragForceMultiplierSinglePoint = g_dragForceMultiplierSinglePoint_Euler;
	}
	else if (jello->integrator[0] == 'R')
	{
		dragForceMultiplierSinglePoint = g_dragForceMultiplierSinglePoint_RK4;
	}
	else if (jello->integrator[0] == 'M')
	{
		dragForceMultiplierSinglePoint = g_dragForceMultiplierSinglePoint_Midpoint;
	}
	struct point sumForce;
	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				//pSUM(fHookLinear[i][j][k], fDampLinear[i][j][k], sumForce);
				pMAKE(
					fHookLinear[i][j][k].x + fDampLinear[i][j][k].x + fCollision[i][j][k].x + fExtForce[i][j][k].x,
					fHookLinear[i][j][k].y + fDampLinear[i][j][k].y + fCollision[i][j][k].y + fExtForce[i][j][k].y,
					fHookLinear[i][j][k].z + fDampLinear[i][j][k].z + fCollision[i][j][k].z + fExtForce[i][j][k].z,
					sumForce);
				if (g_iPickingAMassPoint) // only apply the drag force on this certain point
				{
					if (i == g_pickedPointIndices[0] && j == g_pickedPointIndices[1] && k == g_pickedPointIndices[2])
					{
						struct point mouseDragForceOnSinglePoint;
						pMULTIPLY(g_pMouseDragForce, dragForceMultiplierSinglePoint, mouseDragForceOnSinglePoint);
						pSUM(sumForce, mouseDragForceOnSinglePoint, sumForce);
					}
				}
				else // apply drag force on each point if no point is being selected
				{
					pSUM(sumForce, g_pMouseDragForce, sumForce);
				}
				pMULTIPLY(sumForce, (1.0 / jello->mass), a[i][j][k]);
			}
		}
	}

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

/* performs one step of Euler Midpoint Integration */
/* as a result, updates the jello structure */
void EulerMidpoint(struct world* jello)
{
	int i, j, k;
	point a[8][8][8];

	static struct world buffer;
	static struct point midP[8][8][8], midV[8][8][8];

	buffer = *jello;
	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// fMid
				pMULTIPLY(jello->v[i][j][k], jello->dt, midP[i][j][k]);
				pMULTIPLY(a[i][j][k], jello->dt, midV[i][j][k]);
				pMULTIPLY(midP[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(midV[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	// evaluate f at the mid point
	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// x(t) + deltaT*fmid
				jello->p[i][j][k].x += jello->dt * buffer.v[i][j][k].x;
				jello->p[i][j][k].y += jello->dt * buffer.v[i][j][k].y;
				jello->p[i][j][k].z += jello->dt * buffer.v[i][j][k].z;
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

bool isValidIndex(int i, int j, int k)
{
	return (i >= 0) && (i <= 7) && (j >= 0) && (j <= 7) && (k >= 0) && (k <= 7);
}

void computeLinearHookAndDamp(world* jello, point fHookLinear[8][8][8], point fDampLinear[8][8][8])
{
	// Hook's law: double f_hook=kx
	// rest length: structural 1.0/7, shear: 1.0 / 7 * sqrt(2) OR 1.0 / 7 * sqrt(3), bend: 1.0 / 7 * 2
	static const double restLengthStructural = 1.0 / 7;
	static const double restLengthShearDiagonal2D = restLengthStructural * sqrt(2.0);
	static const double restLengthShearDiagonal3D = restLengthStructural * sqrt(3.0);
	static const double restLengthBend = restLengthStructural * 2;

	// for both linear Hook and linear damp, the forces on the spring are always opposite on the two ends
	// so we need to update both: cache[x1][y1][z1][x2][y2][z2] == -cache[x2][y2][z2][x1][y1][z1] == the value
	static bool fCacheCalculatedHook[8][8][8][8][8][8];
	static bool fCacheCalculatedDamp[8][8][8][8][8][8];
	static struct point fCacheForceHookLinear[8][8][8][8][8][8];
	static struct point fCacheForceDampLinear[8][8][8][8][8][8];

	// Linear Hook: double f_hookLinear = jello->kElastic * x
	// Linear damping: double f_linearDamp=-k_damp * v, where v is the relative speed ON THE SPRING'S DIRECTION
	memset(fCacheCalculatedHook, 0, sizeof(bool) * 512 * 512);
	memset(fCacheCalculatedDamp, 0, sizeof(bool) * 512 * 512);
	memset(fCacheForceHookLinear, 0, sizeof(struct point) * 512 * 512);
	memset(fCacheForceDampLinear, 0, sizeof(struct point) * 512 * 512);

	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				// structural
				{
					if (isValidIndex(i - 1, j, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j, k,
								restLengthStructural, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j, k,
								restLengthStructural, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j - 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j - 1, k,
								restLengthStructural, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j - 1, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j + 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j + 1, k,
								restLengthStructural, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j + 1, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j, k - 1,
								restLengthStructural, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j, k + 1,
								restLengthStructural, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
				}
				// shear2D
				{
					// k-plane
					if (isValidIndex(i - 1, j - 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j - 1, k,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j - 1, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i - 1, j + 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j + 1, k,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j + 1, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j - 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j - 1, k,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j - 1, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j + 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j + 1, k,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j + 1, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					// j-plane
					if (isValidIndex(i - 1, j, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j, k - 1,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i - 1, j, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j, k + 1,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j, k - 1,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j, k + 1,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					// i-plane
					if (isValidIndex(i, j - 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j - 1, k - 1,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j - 1, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j - 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j - 1, k + 1,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j - 1, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j + 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j + 1, k - 1,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j + 1, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j + 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j + 1, k + 1,
								restLengthShearDiagonal2D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j + 1, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
				}
				// shear3D
				{
					// ---
					if (isValidIndex(i - 1, j - 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j - 1, k - 1,
								restLengthShearDiagonal3D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j - 1, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// --+
					if (isValidIndex(i - 1, j - 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j - 1, k + 1,
								restLengthShearDiagonal3D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j - 1, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// -++
					if (isValidIndex(i - 1, j + 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j + 1, k + 1,
								restLengthShearDiagonal3D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j + 1, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// -+-
					if (isValidIndex(i - 1, j + 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j + 1, k - 1,
								restLengthShearDiagonal3D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j + 1, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// +--
					if (isValidIndex(i + 1, j - 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j - 1, k - 1,
								restLengthShearDiagonal3D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j - 1, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// +-+
					if (isValidIndex(i + 1, j - 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j - 1, k + 1,
								restLengthShearDiagonal3D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j - 1, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// +++
					if (isValidIndex(i + 1, j + 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j + 1, k + 1,
								restLengthShearDiagonal3D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j + 1, k + 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// ++-
					if (isValidIndex(i + 1, j + 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j + 1, k - 1,
								restLengthShearDiagonal3D, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j + 1, k - 1,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
				}
				// bend
				{
					if (isValidIndex(i - 2, j, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 2, j, k,
								restLengthBend, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 2, j, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 2, j, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 2, j, k,
								restLengthBend, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 2, j, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j - 2, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j - 2, k,
								restLengthBend, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j - 2, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j + 2, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j + 2, k,
								restLengthBend, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j + 2, k,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j, k - 2)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j, k - 2,
								restLengthBend, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j, k - 2,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j, k + 2)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j, k + 2,
								restLengthBend, fCacheCalculatedHook, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j, k + 2,
								fCacheCalculatedDamp, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
				}
			}
		}
	}
}

struct point computeHookForceOnA(struct world* jello, int xA, int yA, int zA, int xB, int yB, int zB, double restLength,
	bool fCacheForceHookCalculated[8][8][8][8][8][8], struct point fCacheForceHookLinear[8][8][8][8][8][8])
{
	// if (calculated) {fTmp = fCache[i][j][k][i-1][j][k];}
					// else { fTmp = fCache[i][j][k][i-1][j][k] = jello->kElastic kElastic * (x - x0); }
					// fHookLinear[i][j][k] += fTmp;
	struct point fHookResult;
	if (fCacheForceHookCalculated[xA][yA][zA][xB][yB][zB]) {
		fHookResult = fCacheForceHookLinear[xA][yA][zA][xB][yB][zB];
	}
	else
	{
		// calc & write cache
		struct point realSpringVector;
		pDIFFERENCE(jello->p[xA][yA][zA], jello->p[xB][yB][zB], realSpringVector);
		double realLength;
		pLENGTH(realSpringVector, realLength);
		double fHook3D = -jello->kElastic * (realLength - restLength); // -k_hook * (L_length - R)
		struct point realSpringVectorNormalized;
		pCPY(realSpringVector, realSpringVectorNormalized);
		double length; // to use pNORMALIZE
		pNORMALIZE(realSpringVectorNormalized);
		// calc result
		pMAKE(fHook3D * realSpringVectorNormalized.x,
			fHook3D * realSpringVectorNormalized.y,
			fHook3D * realSpringVectorNormalized.z,
			fHookResult);
		// update cache for both [PointA,PointB] & [PointB,PointA]
		pCPY(fHookResult, fCacheForceHookLinear[xA][yA][zA][xB][yB][zB]);
		pMAKE(-fHookResult.x, -fHookResult.y, -fHookResult.z, fCacheForceHookLinear[xB][yB][zB][xA][yA][zA]);
		fCacheForceHookCalculated[xA][yA][zA][xB][yB][zB] = true;
		fCacheForceHookCalculated[xB][yB][zB][xA][yA][zA] = true;
	}
	return fHookResult;
}

struct point computeDampForceOnA(struct world* jello, int xA, int yA, int zA, int xB, int yB, int zB,
	bool fCacheForceDampCalculated[8][8][8][8][8][8], struct point fCacheForceDampLinear[8][8][8][8][8][8])
{
	struct point fDampResult;
	if (fCacheForceDampCalculated[xA][yA][zA][xB][yB][zB]) {
		fDampResult = fCacheForceDampLinear[xA][yA][zA][xB][yB][zB];
	}
	else
	{
		struct point /*realSpringVector,*/ realSpringVectorNormalized; // we only need a normalized vector
		pDIFFERENCE(jello->p[xA][yA][zA], jello->p[xB][yB][zB], realSpringVectorNormalized);
		double length; // to use pNORMALIZE
		pNORMALIZE(realSpringVectorNormalized);

		struct point relativeSpeed;
		pDIFFERENCE(jello->v[xA][yA][zA], jello->v[xB][yB][zB], relativeSpeed);
		double relativeSpeedOnSpring;
		DOTPRODUCTp(relativeSpeed, realSpringVectorNormalized, relativeSpeedOnSpring);
		double fDamp3D = -jello->dElastic * relativeSpeedOnSpring;
		pMAKE(fDamp3D * realSpringVectorNormalized.x,
			fDamp3D * realSpringVectorNormalized.y,
			fDamp3D * realSpringVectorNormalized.z,
			fDampResult);

		pCPY(fDampResult, fCacheForceDampLinear[xA][yA][zA][xB][yB][zB]);
		pMAKE(-fDampResult.x, -fDampResult.y, -fDampResult.z, fCacheForceDampLinear[xB][yB][zB][xA][yA][zA]);
		fCacheForceDampCalculated[xA][yA][zA][xB][yB][zB] = true;
		fCacheForceDampCalculated[xB][yB][zB][xA][yA][zA] = true;
	}
	return fDampResult;
}

void computeCollisionHookAndDampFromBoundingBox(world* jello, point fCollision[8][8][8])
{
	// box's 6 sides: <xmin, >xmax, <ymin, >ymax, <zmin, >zmax
	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				const struct point& curPoint = jello->p[i][j][k];
				const struct point& curPointVelocity = jello->v[i][j][k];
				// here we can easily get collision point because the bounding box's planes are all (x/y/z = +/- 2)
				if (curPoint.x < -2) {
					struct point fCollisionHook = computeCollisionHookForce(curPoint, jello->kCollision, -2, curPoint.y, curPoint.z);
					pSUM(fCollisionHook, fCollision[i][j][k], fCollision[i][j][k]);

					struct point fCollisionDamp = computeCollisionDampForce(curPoint, curPointVelocity, jello->dCollision, -2, curPoint.y, curPoint.z);
					pSUM(fCollisionDamp, fCollision[i][j][k], fCollision[i][j][k]);
				}
				if (curPoint.x > 2) {
					struct point fCollisionHook = computeCollisionHookForce(curPoint, jello->kCollision, 2, curPoint.y, curPoint.z);
					pSUM(fCollisionHook, fCollision[i][j][k], fCollision[i][j][k]);

					struct point fCollisionDamp = computeCollisionDampForce(curPoint, curPointVelocity, jello->dCollision, 2, curPoint.y, curPoint.z);
					pSUM(fCollisionDamp, fCollision[i][j][k], fCollision[i][j][k]);
				}
				if (curPoint.y < -2) {
					struct point fCollisionHook = computeCollisionHookForce(curPoint, jello->kCollision, curPoint.x, -2, curPoint.z);
					pSUM(fCollisionHook, fCollision[i][j][k], fCollision[i][j][k]);

					struct point fCollisionDamp = computeCollisionDampForce(curPoint, curPointVelocity, jello->dCollision, curPoint.x, -2, curPoint.z);
					pSUM(fCollisionDamp, fCollision[i][j][k], fCollision[i][j][k]);
				}
				if (curPoint.y > 2) {
					struct point fCollisionHook = computeCollisionHookForce(curPoint, jello->kCollision, curPoint.x, 2, curPoint.z);
					pSUM(fCollisionHook, fCollision[i][j][k], fCollision[i][j][k]);

					struct point fCollisionDamp = computeCollisionDampForce(curPoint, curPointVelocity, jello->dCollision, curPoint.x, 2, curPoint.z);
					pSUM(fCollisionDamp, fCollision[i][j][k], fCollision[i][j][k]);
				}
				if (curPoint.z < -2) {
					struct point fCollisionHook = computeCollisionHookForce(curPoint, jello->kCollision, curPoint.x, curPoint.y, -2);
					pSUM(fCollisionHook, fCollision[i][j][k], fCollision[i][j][k]);

					struct point fCollisionDamp = computeCollisionDampForce(curPoint, curPointVelocity, jello->dCollision, curPoint.x, curPoint.y, -2);
					pSUM(fCollisionDamp, fCollision[i][j][k], fCollision[i][j][k]);
				}
				if (curPoint.z > 2) {
					struct point fCollisionHook = computeCollisionHookForce(curPoint, jello->kCollision, curPoint.x, curPoint.y, 2);
					pSUM(fCollisionHook, fCollision[i][j][k], fCollision[i][j][k]);

					struct point fCollisionDamp = computeCollisionDampForce(curPoint, curPointVelocity, jello->dCollision, curPoint.x, curPoint.y, 2);
					pSUM(fCollisionDamp, fCollision[i][j][k], fCollision[i][j][k]);
				}
			}
		}
	}
}

struct point computeCollisionHookForce(const point& massPoint, double kCollision,
	double collisionPointX, double collisionPointY, double collisionPointZ)
{
	struct point collisionPoint, collisionSpringVector;
	pMAKE(collisionPointX, collisionPointY, collisionPointZ, collisionPoint);
	pDIFFERENCE(massPoint, collisionPoint, collisionSpringVector);
	double collisionSpringLength;
	pLENGTH(collisionSpringVector, collisionSpringLength);
	double fHookCollision3D = -kCollision * (collisionSpringLength /*- 0, restLength == 0*/);

	struct point fHookCollision, collisionSpringVectorNormalized;
	double length; // to use pNORMALIZE
	pCPY(collisionSpringVector, collisionSpringVectorNormalized);
	pNORMALIZE(collisionSpringVectorNormalized);
	pMAKE(fHookCollision3D * collisionSpringVectorNormalized.x,
		fHookCollision3D * collisionSpringVectorNormalized.y,
		fHookCollision3D * collisionSpringVectorNormalized.z,
		fHookCollision);

	return fHookCollision;
}

point computeCollisionDampForce(const point& massPoint, const point& curPointVelocity, double dCollision, double collisionPointX, double collisionPointY, double collisionPointZ)
{
	struct point collisionPoint, collisionSpringVectorNormalized;
	pMAKE(collisionPointX, collisionPointY, collisionPointZ, collisionPoint);
	pDIFFERENCE(massPoint, collisionPoint, collisionSpringVectorNormalized);
	double length; // to use pNORMALIZE
	pNORMALIZE(collisionSpringVectorNormalized);
	// here vB is always zero because it's a point on the static bounding box
	// vA-vB is always vA so we can simply take vA as the relative speed
	// pDIFFERENCE(curPointVelocity, velocityB, relativeSpeedVector) // we dont need this in this case
	double relativeSpeedOnSpring;
	DOTPRODUCTp(curPointVelocity, collisionSpringVectorNormalized, relativeSpeedOnSpring);
	double fCollisionDamp3D = -dCollision * relativeSpeedOnSpring;
	struct point fCollisionDamp;
	pMAKE(fCollisionDamp3D * collisionSpringVectorNormalized.x,
		fCollisionDamp3D * collisionSpringVectorNormalized.y,
		fCollisionDamp3D * collisionSpringVectorNormalized.z,
		fCollisionDamp);

	return fCollisionDamp;
}

void computeExternalForces(struct world* jello, struct point fExtForce[8][8][8])
{
	double forceFieldGridLength = 4.0 / (jello->resolution - 1);
	double offset = 2.0;
	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				const struct point& curPoint = jello->p[i][j][k];
				if (curPoint.x < -2 || curPoint.y < -2 || curPoint.z < -2 ||
					curPoint.x > 2 || curPoint.y > 2 || curPoint.z > 2) {
					continue;
				}
				double forceFieldGridI = (curPoint.x + offset) / forceFieldGridLength;
				double forceFieldGridJ = (curPoint.y + offset) / forceFieldGridLength;
				double forceFieldGridK = (curPoint.z + offset) / forceFieldGridLength;
				int forceGridIndexI = floor(forceFieldGridI);
				int forceGridIndexJ = floor(forceFieldGridJ);
				int forceGridIndexK = floor(forceFieldGridK);
				double alpha = ((curPoint.x + offset) - (forceGridIndexI * forceFieldGridLength)) / forceFieldGridLength;
				double beta = ((curPoint.y + offset) - (forceGridIndexJ * forceFieldGridLength)) / forceFieldGridLength;
				double gamma = ((curPoint.z + offset) - (forceGridIndexK * forceFieldGridLength)) / forceFieldGridLength;


				if (forceGridIndexI == jello->resolution - 1) { forceGridIndexI--; alpha = 1.0; }
				if (forceGridIndexJ == jello->resolution - 1) { forceGridIndexI--; beta = 1.0; }
				if (forceGridIndexK == jello->resolution - 1) { forceGridIndexI--; gamma = 1.0; }

				//struct point fExtForceOnPoint;
				pMAKE(
					(1 - alpha) * (1 - beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ, forceGridIndexK, jello->resolution)].x +
					(alpha) * (1 - beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ, forceGridIndexK, jello->resolution)].x +
					(1 - alpha) * (beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ + 1, forceGridIndexK, jello->resolution)].x +
					(1 - alpha) * (1 - beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ, forceGridIndexK + 1, jello->resolution)].x +
					(1 - alpha) * (beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ + 1, forceGridIndexK + 1, jello->resolution)].x +
					(alpha) * (beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ + 1, forceGridIndexK, jello->resolution)].x +
					(alpha) * (1 - beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ, forceGridIndexK + 1, jello->resolution)].x +
					(alpha) * (beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ + 1, forceGridIndexK + 1, jello->resolution)].x,

					(1 - alpha) * (1 - beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ, forceGridIndexK, jello->resolution)].y +
					(alpha) * (1 - beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ, forceGridIndexK, jello->resolution)].y +
					(1 - alpha) * (beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ + 1, forceGridIndexK, jello->resolution)].y +
					(1 - alpha) * (1 - beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ, forceGridIndexK + 1, jello->resolution)].y +
					(1 - alpha) * (beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ + 1, forceGridIndexK + 1, jello->resolution)].y +
					(alpha) * (beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ + 1, forceGridIndexK, jello->resolution)].y +
					(alpha) * (1 - beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ, forceGridIndexK + 1, jello->resolution)].y +
					(alpha) * (beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ + 1, forceGridIndexK + 1, jello->resolution)].y,

					(1 - alpha) * (1 - beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ, forceGridIndexK, jello->resolution)].z +
					(alpha) * (1 - beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ, forceGridIndexK, jello->resolution)].z +
					(1 - alpha) * (beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ + 1, forceGridIndexK, jello->resolution)].z +
					(1 - alpha) * (1 - beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ, forceGridIndexK + 1, jello->resolution)].z +
					(1 - alpha) * (beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI, forceGridIndexJ + 1, forceGridIndexK + 1, jello->resolution)].z +
					(alpha) * (beta) * (1 - gamma) * jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ + 1, forceGridIndexK, jello->resolution)].z +
					(alpha) * (1 - beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ, forceGridIndexK + 1, jello->resolution)].z +
					(alpha) * (beta) * (gamma)*jello->forceField[getForceFieldIndex(forceGridIndexI + 1, forceGridIndexJ + 1, forceGridIndexK + 1, jello->resolution)].z,

					fExtForce[i][j][k]
				);

			}
		}
	}
}

void computeCollisionHookAndDampFromIncPlane(struct world* jello, struct point fCollision[8][8][8])
{
	int countPositiveSide = 0, countNegativeSide = 0;
	double fxyz[8][8][8] = { 0.0 };
	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				const struct point& curPoint = jello->p[i][j][k];
				fxyz[i][j][k] = jello->a * curPoint.x + jello->b * curPoint.y + jello->c * curPoint.z + jello->d;
				if (fxyz > 0) { countPositiveSide++; }
				if (fxyz < 0) { countNegativeSide++; }
			}
		}
	}
	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				const struct point& curPoint = jello->p[i][j][k];
				const struct point& curPointVelocity = jello->v[i][j][k];
				if (countPositiveSide > countNegativeSide) {
					// should calc all points with negative fxyz
					if (fxyz[i][j][k] >= 0) { continue; } // ignore "positive" points
				}
				else {
					// should calc all points with positive fxyz
					if (fxyz[i][j][k] <= 0) { continue; } // ignore "negative" points
				}

				// collision point: P_projection = x/y/z + phi * a/b/c,
				// where phi = -(ax+by+cz+d)/(a*a+b*b+c*c), (x,y,z) are curPoint's position
				double phi = -fxyz[i][j][k] / (jello->a * jello->a + jello->b * jello->b + jello->c * jello->c);
				//struct point collisionPoint;
				double collisionPointX = curPoint.x + phi * jello->a,
					collisionPointY = curPoint.y + phi * jello->b,
					collisionPointZ = curPoint.z + phi * jello->c;
				struct point fCollisionIncPlaneHookForce = computeCollisionHookForce(curPoint, jello->kCollision, collisionPointX, collisionPointY, collisionPointZ);
				struct point fCollisionIncPlaneDampForce = computeCollisionDampForce(curPoint, curPointVelocity, jello->dCollision, collisionPointX, collisionPointY, collisionPointZ);
				struct point fCollisionForceIncPlane;
				pSUM(fCollisionIncPlaneHookForce, fCollisionIncPlaneDampForce, fCollisionForceIncPlane);
				pSUM(fCollisionForceIncPlane, fCollision[i][j][k], fCollision[i][j][k]);
			}
		}
	}
}
