<Please submit this file with your solution.>

CSCI 520, Assignment 1

Marshall Wang

================

- NOTICE:
  - If you want to directly run the executable and capture the screen shots, please create a directory named as "screenShots", since I modified the file saving paths.
    - You don't need to do that manually if you build & run it at Visual Studio. I added some post-build events taking care of this.
  - If you are to run createWorld.exe, please create a directory named as "world" at the same folder.

<Description of what you have accomplished>
*Note: tested successfully on VS2019*
Everything in the requirements:
- finished computeAcceleration() and can animate the jello cube based on a realistic physical model
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

- implemented mouse-drag force
  - you can drag the mouse with the left button and apply forces on the camera's right & up axes, according to the mouse movement on screen.
  - what I did was to get the camera's viewMatrix by glGetFloatv(GL_MODELVIEW_MATRIX, Glfloat[16])
    - and then extract the first & second column as the camera's Right & Up axes
    - then we can use like mouseDelta to generate the forces and then apply them in computeAcceleration()
  - you can change g_mouseDragForceMultiplier if you have a hard time noticing this effect.

- implemented midpoint method 
  - can be invoked by setting the first line of .w files to "Midpoint" (jello.integrator[0]=='M')
  - I'm not sure if it should be that name but anyway I named the function as "EulerMidpoint"
  - more stable than Euler, but less stable than RK4 (of course)
    - my test: using default configurations of the provided rotate.w, except I changed the timestep
      - when timestep=0.001, only RK4 works
        - Midpoint will not blow up here, but the cube keeps vibrating, which is not correct
      - when timestep=0.0008, both RK4 & Midpoint work, Euler blows up
      - when timestep=0.0005, all the three methods work
