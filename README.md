VennEuler.jl
============

Generate area-proportional Venn/Euler diagrams in Julia. This is based, in part, on 
algorithms from [Leland Wilkinson](http://www.cs.uic.edu/~wilkinson/Publications/venneuler.pdf).

Wilkinson, Leland. "Exact and approximate area-proportional circular Venn and Euler diagrams." _IEEE Transactions on Visualization and Computer Graphics_ 18, no. 2 (2012): 321-331.

GitHub: [HarlanH/VennEuler.jl](https://github.com/HarlanH/VennEuler.jl)

See LICENSE.md for the MIT license.

What It Does
============

Area-proportional Venn/Euler Diagrams show the overlap between sets in such a way that the size of the shapes
is proportional to the size of the sets, and the size of the overlaps on the page
is proportional to the size of the overlaps of the sets. In general, using circles, you can only do
this perfectly if you have two sets -- there will always be some residual error, where the sizes are
not perfectly proportional.

Wilkinson developed a straightforward method of approximately fitting area-proportional diagrams, but
the code was written in Java and was difficult to extend. This Julia package re-implements the 
algorithm, with the following additions:

1. You can use other shapes -- currently squares, triangles, and rectangles, in addition to circles.
2. You can use 3-parameter (X, Y, Q) shapes, such as axis-parallel rectangles, in addition to
   2-parameter (X, Y) shapes, such as squares.
3. There is a relatively easy-to-use specification structure that lets you mix and match shapes.
4. You can lock any parameters you'd like and prevent them from being improved, which is handy for
   putting the largest shape in the middle.
5. It should be easy/easier for others to collaborate and extend.

How It Works
============

You define a structure that specifies the shapes and additional constraints you'd like on the fitting. That
gets translated into a "state" vector, which has 2 (or more) values per shape. A set of lower and upper 
bounds is defined that keeps the shapes all inside the lines. Then the state vector and bounds are 
optimized, currently with a rather brute-force optimizer. To compute the cost function, the shapes are 
drawn onto in-memory bitmaps, the number of overlapping pixels is counted, and the distance between
the resulting counts and the actual target overlap vector is computed. That cost function is minimized.

The results can be rendered as an SVG file.

How To Use It
=============

Dependencies:

* IterTools.jl -- handy!
* NLopt.jl -- the global optimization library (I don't love it, but it seems to work)
* Cairo -- graphics, for output

See the `test/test.jl` and `test/DC2.jl` scripts for examples. The latter shows (as of Spring 2014)
the overlap in membership of six [Data Community DC](http://datacommunitydc.org) Meetup groups.

How To Make It Better
=====================

Pull requests welcome! See the [issues list](https://github.com/HarlanH/VennEuler.jl/issues?state=open)!

