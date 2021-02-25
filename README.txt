<Please submit this file with your solution.>

CSCI 520, Assignment 1

Marshall Wang

================
*Note: test was successful on VS2019, I did not test the project files under other IDE/platforms*
*You should be able to build the solution without any further concerns since I fixed some project settings"

- NOTICE:
  - for executables:
    - If you want to directly run the executable and capture the screen shots, please create a directory named as "screenShots", since I modified the file saving paths.
      - You don't need to do that manually if you build & run it at Visual Studio. I added some post-build events taking care of this.
    - If you are to run createWorld.exe, please create a directory named as "world" at the same folder. It will create a my_jello.w file there
  - world files: 
    - I modified the rotate.w since it always blows up when using Euler. The new file being rotate_my.w, can be simulated stably by Euler.
      - In rotate_my.w, I changed both timestep and the collision/damping coefficients to make sure it's always stable in most cases
      - I also implemented Midpoint method, you can specify 'Midpoint' at the first line of world files to invoke it. Please see the last section for details.
      - rotate_my.w should always work in all the three integrator mode. (Euler/Midpoint/RK4)

<Description of what you have accomplished>
I finished everything in the core requirements:
- computeAcceleration() and can animate the jello cube based on a realistic physical model
- collision response
  - with bounding box
    - since the bounding box is stationed there, I can simply compare each point's x/y/z-axis positions to judge if it's colliding with the six sides
    - and it's fairly easy to generate collision spring because the collision point is also easy to find & calculate
  - with inclined plane (see the last section)
- run at above 15fps under both Euler & RK4 when n = 5. 
  - minimized memory allocation
  - cached the results when calculating linear hook & damping
  - CPU: i7-9900K, 20% usage when running

Misc:
- changed some VS project settings & added some post-build events, now they can be compiled and built without any errors
- changed lighting/material settings to my preferred color theme, cyan-ish
- tested simple animation: move each mass points by its (jello.v * dt)
- restructured doIdle(), where inside it calls Animate(), captureScreenShots() and then glutPostRedisplay()
- changed macro pMAKE, making it more compatible with struct point
- changed the description of pDIFFERENCE, which should be (src1 - src2)
- added macro pLENGTH

// -----------------------------------------------------------------------------------------
<Also, explain any extra credit that you have implemented.>
- collision with inclined plane
  - collision point: it's basically solving the equation: F(P_outplane + alpha * vecPlaneNormal) = 0 for alpha
  - which part of the cube is colliding:
    - what I was doing was to count the points lying on each side of the plane
    - then the minority of points should be our results

- implemented midpoint method 
  - can be invoked by setting the first line of .w files to "Midpoint" (jello.integrator[0]=='M')
  - I'm not sure if it should be that name but anyway I named the function as "EulerMidpoint"
  - more stable than Euler, but less stable than RK4 (of course)
    - my test: using default configurations of the provided rotate.w, except I changed the timestep
      - when timestep=0.001, only RK4 works
        - Midpoint will not blow up here, but the cube keeps vibrating, which is not correct
      - when timestep=0.0008, both RK4 & Midpoint work, Euler blows up
      - when timestep=0.0005, all the three methods work

- implemented mouse-drag force
  - this includes both applying forces to every point (when you are not clicking on a mass point) and to some single points (when you are dragging a mass point)
  - you can drag the mouse with the left button and apply forces on the camera's right & up axes, according to the mouse movement on screen.
  - How I generate the force:
    - what I did was to get the camera's viewMatrix by glGetFloatv(GL_MODELVIEW_MATRIX, Glfloat[16])
      - and then extract the first & second column as the camera's Right & Up axes
      - then we can use like mouseDelta to generate the forces and then apply them in computeAcceleration()
  - Picking a specific point:
    - It's much easier to select a point when in polygon mode than in wireframe mode
      - I used a stencil buffer, pixels that contains the cube will have stencilValue == 1
      - Then what I did was basically get the mouseClick point and called gluUnproject() to transform it back to world space
      - so that I can easily find and attach a closest point
    - NOTE: If you are dragging a single point, there will be no mouse-dragging forces being applied to any other points
  - I used different force magnitude (actually multiplier) in each different mode:
    - RK4 allows the largest magnitude, then Midpoint, and Euler has the worse case
    - If you are simply applying forces to every point by dragging in the blanks, the forces is not multiplied
  - You can change these settings in jello.cpp between line 31 and line 35.

