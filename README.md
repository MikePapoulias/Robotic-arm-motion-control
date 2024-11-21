Motion control of a UR10e robotic arm with 6 degrees of freedom and sphere tracking on a curved surface.

# Description

In this project, we address a problem involving a UR10e-6 robotic arm and a curved surface on which a sphere moves. The task is to design an appropriate control signal for the joint velocities of the arm, ensuring that its end-effector follows the moving sphere while maintaining a specific relative position and orientation. The problem's simulation and the required graphs are implemented in the Matlab environment.

The robotic arm starts with predefined initial joint angles and is subject to specific maximum absolute velocity and acceleration limits for each joint. The sphere moves under the influence of gravity on the curved surface and decelerates continuously due to friction.

This process is tackled in Part A, where we design the control signal to track the sphere. In Part B, the objective shifts to designing a control signal for the joint velocities so that the robot's end-effector, now equipped with a gripper, successfully grasp the moving sphere.

Required time-dependent plots include, among others, the position, velocity, and acceleration responses for each joint individually.
