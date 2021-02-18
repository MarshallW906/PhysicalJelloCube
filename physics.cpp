/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

static bool isValidIndex(int i, int j, int k);
static struct point computeHookForceOnA(struct world* jello, int xA, int yA, int zA, int xB, int yB, int zB, double restLength,
	bool fCacheForceCalculated[8][8][8][8][8][8], struct point fCacheForceHookLinear[8][8][8][8][8][8]);
static struct point computeDampForceOnA(struct world* jello, int xA, int yA, int zA, int xB, int yB, int zB,
	bool fCacheForceCalculated[8][8][8][8][8][8], struct point fCacheForceDampLinear[8][8][8][8][8][8]);

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
	static struct point fCacheForceHookLinear[8][8][8][8][8][8];
	static struct point fCacheForceDampLinear[8][8][8][8][8][8];

	// Linear Hook: double f_hookLinear = jello->kElastic * x
	// Linear damping: double f_linearDamp=-k_damp * v, where v is the relative speed ON THE SPRING'S DIRECTION
	memset(fCacheForceCalculated, 0, sizeof(bool) * 512 * 512);
	memset(fCacheForceHookLinear, 0, sizeof(struct point) * 512 * 512);
	memset(fCacheForceDampLinear, 0, sizeof(struct point) * 512 * 512);

	struct point fHookLinear[8][8][8] = { {{0.0}} };
	struct point fDampLinear[8][8][8] = { {{0.0}} };
	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				// structural
				{
					if (isValidIndex(i - 1, j, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j - 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j - 1, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j - 1, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j + 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j + 1, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j + 1, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
				}
				// shear2D
				{
					// k-plane
					if (isValidIndex(i - 1, j - 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j - 1, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j - 1, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i - 1, j + 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j + 1, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j + 1, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j - 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j - 1, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j - 1, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j + 1, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j + 1, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j + 1, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					// j-plane
					if (isValidIndex(i - 1, j, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i - 1, j, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 1, j, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					// i-plane
					if (isValidIndex(i, j - 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j - 1, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j - 1, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j - 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j - 1, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j - 1, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j + 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j + 1, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j + 1, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j + 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j + 1, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j + 1, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
				}
				// shear3D
				{
					// ---
					if (isValidIndex(i - 1, j - 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j - 1, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j - 1, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// --+
					if (isValidIndex(i - 1, j - 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j - 1, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j - 1, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// -++
					if (isValidIndex(i - 1, j + 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j + 1, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j + 1, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// -+-
					if (isValidIndex(i - 1, j + 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 1, j + 1, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 1, j + 1, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// +--
					if (isValidIndex(i + 1, j - 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j - 1, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j - 1, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// +-+
					if (isValidIndex(i + 1, j - 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j - 1, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j - 1, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// +++
					if (isValidIndex(i + 1, j + 1, k + 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j + 1, k + 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j + 1, k + 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
					// ++-
					if (isValidIndex(i + 1, j + 1, k - 1)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 1, j + 1, k - 1,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 1, j + 1, k - 1,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
				}
				// bend
				{
					if (isValidIndex(i - 2, j, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i - 2, j, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i - 2, j, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i + 2, j, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i + 2, j, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i + 2, j, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j - 2, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j - 2, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j - 2, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j + 2, k)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j + 2, k,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j + 2, k,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j, k - 2)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j, k - 2,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j, k - 2,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}

					if (isValidIndex(i, j, k + 2)) {
						struct point fHookResult =
							computeHookForceOnA(jello, i, j, k, i, j, k + 2,
								restLengthStructural, fCacheForceCalculated, fCacheForceHookLinear);
						pSUM(fHookLinear[i][j][k], fHookResult, fHookLinear[i][j][k]);
						struct point fDampResult =
							computeDampForceOnA(jello, i, j, k, i, j, k + 2,
								fCacheForceCalculated, fCacheForceDampLinear);
						pSUM(fDampLinear[i][j][k], fDampResult, fDampLinear[i][j][k]);
					}
				}
			}
		}
	}


	// Force field: double f_extForceField
	// trilinear interpolate and get the actual force at certain points
	//struct point fExtForce[8][8][8];

	// collision

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

bool isValidIndex(int i, int j, int k)
{
	return (i >= 0) && (i <= 7) && (j >= 0) && (j <= 7) && (k >= 0) && (k <= 7);
}

struct point computeHookForceOnA(struct world* jello, int xA, int yA, int zA, int xB, int yB, int zB, double restLength,
	bool fCacheForceCalculated[8][8][8][8][8][8], struct point fCacheForceHookLinear[8][8][8][8][8][8])
{
	// if (calculated) {fTmp = fCache[i][j][k][i-1][j][k];}
					// else { fTmp = fCache[i][j][k][i-1][j][k] = jello->kElastic kElastic * (x - x0); }
					// fHookLinear[i][j][k] += fTmp;
	struct point fHookResult;
	if (fCacheForceCalculated[xA][yA][zA][xB][yB][zB]) {
		fHookResult = fCacheForceHookLinear[xA][yA][zA][xB][yB][zB];
	}
	else {
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
		fCacheForceCalculated[xA][yA][zA][xB][yB][zB] = true;
		fCacheForceCalculated[xB][yB][zB][xA][yA][zA] = true;
	}
	return fHookResult;
}

struct point computeDampForceOnA(world* jello, int xA, int yA, int zA, int xB, int yB, int zB,
	bool fCacheForceCalculated[8][8][8][8][8][8], point fCacheForceDampLinear[8][8][8][8][8][8])
{
	// first project both vA & vB onto the target spring by dotProduct(v, dirNormalized) 
	// then calc
	return point();
}
