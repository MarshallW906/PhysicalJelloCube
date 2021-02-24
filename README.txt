<Please submit this file with your solution.>

CSCI 520, Assignment 1

Marshall Wang

================

<Description of what you have accomplished>
*Note: tested on VS2019, not yet tested for VS2017*
- collision with inclined plane
- implemented midpoint method 
  - I'm not sure if it should be that name but anyway I named the function as "EulerMidpoint"
  - can be invoked by setting the first line of .w files to "Midpoint" (jello.integrator[0]=='M')
- implemented mouse-drag force
 - you can drag the mouse with the left button and apply forces on the camera's right & up axes, according to the mouse movement on screen.

*My personal progress tracker*
What I have done:
- changed some VS project settings & added some post-build events, now they can be compiled and built without any errors
- fixed the position error at createWorld.cpp
- changed some lighting/material settings to my preferred color, cyan-ish
- tested simple animation: move each mass points by its (jello.v * dt)
- restructured doIdle(), where inside it calls Animate(), captureScreenShots() and then glutPostRedisplay()
- changed macro pMAKE, making it more compatible with struct point


<Also, explain any extra credit that you have implemented.>

