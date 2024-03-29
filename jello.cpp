/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

  Your name:
  Marshall Wang

*/

#include "jello.h"
#include "showCube.h"
#include "input.h"
#include "physics.h"

#include <chrono>   

// camera parameters
double Theta = pi / 6;
double Phi = pi / 6;
double R = 6;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton, g_iMiddleMouseButton, g_iRightMouseButton;

// number of images saved to disk so far
int sprite = 0;

double g_mouseDragForceBaseMultiplier = 0.02;
struct point g_pMouseDragForce = { 0.0 };
double g_dragForceMultiplierSinglePoint_Euler = 10;
double g_dragForceMultiplierSinglePoint_Midpoint = 40;
double g_dragForceMultiplierSinglePoint_RK4 = 400;
int g_iPickingAMassPoint = 0;
int g_pickedPointIndices[3] = { -1, -1, -1 };

// these variables control what is displayed on screen
int shear = 0, bend = 0, structural = 1, pause = 0, viewingMode = 0, saveScreenToFile = 0;

struct world jello;

int windowWidth, windowHeight;

std::chrono::system_clock::time_point lastFrameTimePoint;
double deltaTime = 0.0;

void myinit()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0, 1.0, 0.01, 1000.0);

	// set background color to grey
	glClearColor(0.5, 0.5, 0.5, 0.0);
	glClearStencil(0);

	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);

	glShadeModel(GL_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	lastFrameTimePoint = std::chrono::system_clock::now();

	return;
}

void reshape(int w, int h)
{
	// Prevent a divide by zero, when h is zero.
	// You can't make a window of zero height.
	if (h == 0)
		h = 1;

	glViewport(0, 0, w, h);

	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// Set the perspective
	double aspectRatio = 1.0 * w / h;
	gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	windowWidth = w;
	windowHeight = h;

	glutPostRedisplay();
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// camera parameters are Phi, Theta, R
	gluLookAt(R * cos(Phi) * cos(Theta), R * sin(Phi) * cos(Theta), R * sin(Theta),
		0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

	/* Lighting */
	/* You are encouraged to change lighting parameters or make improvements/modifications
	   to the lighting model .
	   This way, you will personalize your assignment and your assignment will stick out.
	*/

	// global ambient light
	// now ambient: none (black)
	GLfloat aGa[] = { 0.0, 0.0, 0.0, 0.0 };

	// light 's ambient, diffuse, specular
	GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd0[] = { 1.0, 1.0, 1.0, 1.0 }; // white
	GLfloat lKs0[] = { 1.0, 1.0, 1.0, 1.0 }; // white

	GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd1[] = { 1.0, 0.0, 0.0, 1.0 }; // red
	GLfloat lKs1[] = { 1.0, 0.0, 0.0, 1.0 }; // red

	GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd2[] = { 1.0, 1.0, 0.0, 1.0 }; // yellow
	GLfloat lKs2[] = { 1.0, 1.0, 0.0, 1.0 }; // yellow

	GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd3[] = { 0.0, 1.0, 1.0, 1.0 }; // cyan
	GLfloat lKs3[] = { 0.0, 1.0, 1.0, 1.0 }; // cyan

	GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd4[] = { 0.0, 0.0, 1.0, 1.0 }; // blue
	GLfloat lKs4[] = { 0.0, 0.0, 1.0, 1.0 }; // blue

	GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd5[] = { 1.0, 0.0, 1.0, 1.0 }; // purple/magenta
	GLfloat lKs5[] = { 1.0, 0.0, 1.0, 1.0 }; // purple/magenta

	GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd6[] = { 1.0, 1.0, 1.0, 1.0 }; // white
	GLfloat lKs6[] = { 1.0, 1.0, 1.0, 1.0 }; // white

	GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd7[] = { 0.0, 1.0, 1.0, 1.0 }; // cyan
	GLfloat lKs7[] = { 0.0, 1.0, 1.0, 1.0 }; // cyan

	// light positions and directions
	// so by default they are 8 point lights stationed at the 8 sections divided by the xyz axes
	// if w==0, then the light is treated as a directional light and attenuation is disabled
	// otherwise it's a point light
	// lemme change some of them to directional lights
	/* default values
	GLfloat lP0[] = { -1.999, -1.999, -1.999, 1.0 };
	GLfloat lP1[] = { 1.999, -1.999, -1.999, 1.0 };
	GLfloat lP2[] = { 1.999, 1.999, -1.999, 1.0 };
	GLfloat lP3[] = { -1.999, 1.999, -1.999, 1.0 };
	GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
	GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
	GLfloat lP6[] = { 1.999, 1.999, 1.999, 1.0 };
	GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };
	*/
	GLfloat lP0[] = { -1.999, -1.999, -1.999, 0.0 };
	GLfloat lP1[] = { 1.999, -1.999, -1.999, 0.0 };
	GLfloat lP2[] = { 1.999, 1.999, -1.999, 0.0 };
	GLfloat lP3[] = { -1.999, 1.999, -1.999, 0.0 };
	GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
	GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
	GLfloat lP6[] = { 1.999, 1.999, 1.999, 1.0 };
	GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };

	// jelly material color

	GLfloat mKa[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat mKd[] = { 0.3, 0.3, 0.3, 1.0 }; // provided default value: 0.7*black + 0.3*white "grey"
	//GLfloat mKd[] = { 0.2, 0.7, 0.7, 1.0 };   // lets make it a cyan-ish cube
	//GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 }; // default: white
	GLfloat mKs[] = { 1.0, 0.0, 0.0, 1.0 }; // a red specular color 
	//GLfloat mKe[] = { 0.0, 0.0, 0.0, 1.0 }; // emission, default: none(black)
	GLfloat mKe[] = { 0.0, 0.25, 0.25, 1.0 }; // cyan-ish 

	/* set up lighting */
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

	// set up cube color
	glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
	glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
	glMaterialf(GL_FRONT, GL_SHININESS, 10); // default is 120

	// macro to set up light i
#define LIGHTSETUP(i)\
  glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
  glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
  glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
  glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
  glEnable(GL_LIGHT##i)

	// it CAN support more than eight lights
	// but we will need to use "GL_LIGHT0 + n", where n is the n-th light we are setting
	// 8 lights per polygon is plenty though
	LIGHTSETUP(0);
	LIGHTSETUP(1);
	LIGHTSETUP(2);
	LIGHTSETUP(3);
	LIGHTSETUP(4);
	LIGHTSETUP(5);
	LIGHTSETUP(6);
	LIGHTSETUP(7);

	// enable lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_STENCIL_TEST); // Enables testing AND writing functionalities
	//glStencilMask(0xFF); // 0xFF is the default value

	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

	// show the cube
	// inside it writes 1 to the stencil when rendering the jello cube, both wireframe & polygon mode
	showCube(&jello);

	// disable stencil buffer modifications
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

	glDisable(GL_LIGHTING);

	// show the bounding box
	showBoundingBox();

	// show inclined plane
	showIncPlaneIfExists(&jello);

	glutSwapBuffers();
}

void doIdle();
void animateWithoutPhysics();
void animateByPhysics();
void printFrameRate();
void captureScreenShots();


int main(int argc, char** argv)
{
	if (argc < 2)
	{
		printf("Oops! You didn't say the jello world file!\n");
		printf("Usage: %s [worldfile]\n", argv[0]);
		exit(0);
	}

	readWorld(argv[1], &jello);

	glutInit(&argc, argv);

	/* double buffered window, use depth testing, 640x480 */
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	windowWidth = 640;
	windowHeight = 480;
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Jello cube");

	/* tells glut to use a particular display function to redraw */
	glutDisplayFunc(display);

	/* replace with any animate code */
	// I re-structured doIdle() and now it calls physics/simple animation from inside
	glutIdleFunc(doIdle);

	/* callback for mouse drags */
	glutMotionFunc(mouseMotionDrag);

	/* callback for window size changes */
	glutReshapeFunc(reshape);

	/* callback for mouse movement */
	glutPassiveMotionFunc(mouseMotion);

	/* callback for mouse button changes */
	glutMouseFunc(mouseButton);

	/* register for keyboard events */
	glutKeyboardFunc(keyboardFunc);

	/* do initialization */
	myinit();

	/* forever sink in the black hole */
	glutMainLoop();

	return(0);
}

void doIdle()
{
	printFrameRate();

	captureScreenShots();

	if (pause == 0) 
	{
		//animateWithoutPhysics();
		animateByPhysics();
	}

	glutPostRedisplay();
}

void animateWithoutPhysics()
{
	static const double speed = 2.5;

	for (int i = 0; i <= 7; i++) {
		for (int j = 0; j <= 7; j++) {
			for (int k = 0; k <= 7; k++) {
				jello.p[i][j][k].z += jello.dt * jello.v[i][j][k].z * speed;
				jello.p[i][j][k].y += jello.dt * jello.v[i][j][k].y * speed;
				jello.p[i][j][k].x += jello.dt * jello.v[i][j][k].x * speed;
			}
		}
	}
}

void animateByPhysics()
{
	for (int i = 1; i <= jello.n; i++)
	{
		if (jello.integrator[0] == 'E') // Euler
		{
			Euler(&jello);
		}
		if (jello.integrator[0] == 'R') // RK4
		{
			RK4(&jello);
		}
		if (jello.integrator[0] == 'M') // Midpoint
		{
			EulerMidpoint(&jello);
		}
	}
}

void printFrameRate()
{
	auto nowTime = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(nowTime - lastFrameTimePoint);
	deltaTime = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
	double frameRate = 1.0 / deltaTime;
	lastFrameTimePoint = nowTime;
	printf("frame rate: %.2lf\n", frameRate);
}

void captureScreenShots()
{
	char s[50] = "screenShots\\picxxxx.ppm";
	int i;

	// save screen to file
	/*s[3] = 48 + (sprite / 1000);
	s[4] = 48 + (sprite % 1000) / 100;
	s[5] = 48 + (sprite % 100) / 10;
	s[6] = 48 + sprite % 10;*/

	s[15] = 48 + (sprite / 1000);
	s[16] = 48 + (sprite % 1000) / 100;
	s[17] = 48 + (sprite % 100) / 10;
	s[18] = 48 + sprite % 10;

	if (saveScreenToFile == 1)
	{
		saveScreenshot(windowWidth, windowHeight, s);
		//saveScreenToFile = 0; // save only once, change this if you want continous image generation (i.e. animation)
		sprite++;
	}

	if (sprite >= 300) // allow only 300 snapshots
	{
		exit(0);
	}

	/*
	// I moved it into doIdle()
	if (pause == 0)
	{
		// insert code which appropriately performs one step of the cube simulation:
	}
	*/
}