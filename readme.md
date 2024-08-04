LIB INERTIAL
------------

Inertial navigation library functions written in MATLAB, C, and C++ for
research and embedded projects.

Author: Jan Zwiener (jan@zwiener.org)


Features
--------

 - Developed with embedded applications in mind. Worst-case execution time, stack usage.
   Avoid dynamic memory allocation.


C++ Version with Eigen Math Library
-----------------------------------

The Eigen C++ math library is fantastic and the math library that I enjoy the most to work with. No complex dependency, very stable, you can just copy the Eigen math library in your repository if you want to. The resulting code dealing with matrices is quite readable and probably the best (and maybe only valid) use of C++ operator overloading in any C++ project that I've ever seen. The fascinating aspect for me is that the generated code is super fast as well.
That being said, I still would not use it for safety critical embedded applications. The C++ implementation of Eigen is hard to read for me, to be honest, I don't fully understand it, which is not the best foundation for a safety critical application. In those embedded safety critical applications the clear winner for me is the C implementation of the UDU ("Bierman/Thornton") Kalman filter.



