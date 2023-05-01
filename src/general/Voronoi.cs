using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Godot.Collections;
using Ray = MathUtils.Ray;
using Vector3 = Godot.Vector3;

public unsafe class Voronoi
{
    public ICollection<Cell>? VoronoiDiagram;
    public List<Tetrahedron> DelaunayDiagram;
    public Array<Triangle>? Mesh;

    // very big tetrahedron, will automate making this
    private readonly float[][] bigTetra =
    {
        new[] { 0, -24999.75f, -74999.25f },
        new[] { 0, 74999.25f, 0 },
        new[] { 53032.61968f, -24999.75f, -53032.61968f },
        new[] { -53032.61968f, -24999.75f, -53032.61968f },
    };

    public Voronoi(List<Vector3> seeds)
    {
        // big tetrahedron that contains all points to start
        var bigTetVerts = new[]
        {
            new Vector3(bigTetra[0][0], bigTetra[0][1], bigTetra[0][2]),
            new Vector3(bigTetra[1][0], bigTetra[1][1], bigTetra[1][2]),
            new Vector3(bigTetra[2][0], bigTetra[2][1], bigTetra[2][2]),
            new Vector3(bigTetra[3][0], bigTetra[3][1], bigTetra[3][2]),
        };

        var bigTet = new Tetrahedron(bigTetVerts);
        DelaunayDiagram = new List<Tetrahedron> { bigTet };

        InitializeDiagram(seeds);
        DelIso(DelaunayDiagram);
    }

    /// <summary>
    ///   <para>
    ///     Hugo Ledoux 'Computing the 3D Voronoi Diagram Robustly: An Easy Explanation'
    ///   </para>
    ///   Delft University of Technology (OTB-section GIS Technology) [Internet] 2007
    ///   <para>
    ///     Available from: http://www.gdmc.nl/publications/2007/Computing_3D_Voronoi_Diagram.pdf
    ///   </para>
    /// </summary>
    private void InitializeDiagram(List<Vector3> seeds)
    {
        // insert seeds as query points, then rebuild diagram until we triangulate every point.
        // gives low poly initial triangulation that we make more detailed with Del-Iso
        for (int i = 0; i < seeds.Count; i++)
        {
            Vector3 query = seeds[i];
            InsertPoint(DelaunayDiagram, query);
        }
    }

    // this inserts a point and calculates new tetrahedra
    private void InsertPoint(List<Tetrahedron> delaunay, Vector3 point)
    {
        var lastTetra = delaunay.Count - 1;
        Tetrahedron tetra = Walk(delaunay[lastTetra], point);

        DelaunayDiagram.Remove(delaunay[lastTetra]);

        // insert point in tetra with a flip14
        Tetrahedron a = new(new[] { point, tetra.Vertices[0], tetra.Vertices[1], tetra.Vertices[2] });
        Tetrahedron b = new(new[] { point, tetra.Vertices[0], tetra.Vertices[2], tetra.Vertices[3] });
        Tetrahedron c = new(new[] { point, tetra.Vertices[0], tetra.Vertices[3], tetra.Vertices[1] });
        Tetrahedron d = new(new[] { point, tetra.Vertices[1], tetra.Vertices[2], tetra.Vertices[3] });

        Queue<Tetrahedron> newTetras = new();
        newTetras.Enqueue(a);
        newTetras.Enqueue(b);
        newTetras.Enqueue(c);
        newTetras.Enqueue(d);

        while (newTetras.Count > 0)
        {
            // tetra = {p, a, b, c} <--pop from stack
            var test = newTetras.Dequeue();

            // tetra[a] = {a, b, c, d} <--get adjacent tetra of delaunay having abc as a face
            var neighbor = *test.Neighbors[3];

            // if d is inside circumsphere of tetra then flip
            if (InSphere(test.Vertices[0], test.Vertices[1], test.Vertices[2], test.Vertices[3],
                    neighbor.Vertices[0]) > 0)
            {
                Flip(test, neighbor, ref newTetras);
            }
        }

        DelaunayDiagram.Add(a);
        DelaunayDiagram.Add(b);
        DelaunayDiagram.Add(c);
        DelaunayDiagram.Add(d);
    }

    /// <summary>
    ///   <para>
    ///     O.Devillers, S.Pion, and M.Teillaud 'Walking in a Triangulation'
    ///   </para>
    ///   International Journal of Foundations of Computer Science, 13(2):181-199, 2002
    ///   <para>
    ///     Available from: https://inria.hal.science/inria-00102194/document
    ///   </para>
    /// </summary>
    /// <param name="tetra">tetrahedron to walk from</param>
    /// <param name="point">query point</param>
    private Tetrahedron Walk(Tetrahedron tetra, Vector3 point)
    {
        // Remembering Stochastic Walk
        // from some starting vertex to point. tetra=qrlt is a tetrahedron of simplexes (faces)
        // previous = tetra; end = false;
        var previous = tetra;
        bool end = false;

        while (!end)
        {
            // face = random facet of tetra;
            int random = new Random().Next(3);
            Triangle face = tetra.Faces[random];

            // TODO: need to map neighbors to edges
            Tetrahedron neighbor = *tetra.Neighbors[random];

            // don't go backwards
            bool neighborPrevious = neighbor == previous;

            // point on other side of edge
            bool otherSide = Orient(face.Vertices[0], face.Vertices[1], face.Vertices[2], point) > 0;

            // if( point not neighbor of previous through face ) && ( point on other side of face )
            if (!neighborPrevious && otherSide)
            {
                previous = tetra;
                tetra = neighbor;
            }
            else
            {
                // face = next facet of tetra;
                face = random < 3 ? tetra.Faces[random++] : tetra.Faces[0];

                neighbor = *tetra.Neighbors[random];
                neighborPrevious = neighbor == previous;
                otherSide = Orient(face.Vertices[0], face.Vertices[1], face.Vertices[2], point) > 0;

                if (!neighborPrevious && otherSide)
                {
                    previous = tetra;
                    tetra = neighbor;
                }
                else
                {
                    // face = next facet of tetra;
                    face = random < 3 ? tetra.Faces[random++] : tetra.Faces[0];

                    neighbor = *tetra.Neighbors[random];
                    neighborPrevious = neighbor == previous;
                    otherSide = Orient(face.Vertices[0], face.Vertices[1], face.Vertices[2], point) > 0;

                    if (!neighborPrevious && otherSide)
                    {
                        previous = tetra;
                        tetra = neighbor;
                    }
                    else
                    {
                        end = true;
                    }
                }
            }
        }

        // tetra contains point
        return tetra;
    }

    // A technique where we split up a tetrahedron
    // TODO: should map adjacencies here
    // TODO: try testing all edges for relation to edge pd
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void Flip(Tetrahedron tetraA, Tetrahedron tetraB, ref Queue<Tetrahedron> newTetras)
    {
        Ray intersectCheck = new(tetraA.Vertices[0], tetraB.Vertices[0] - tetraA.Vertices[0]);
        Tetrahedron pdab = new(new[]

        {
            tetraA.Vertices[0], tetraB.Vertices[0], tetraA.Vertices[1], tetraA.Vertices[2],
        });

        // case #1: convex
        if (MathUtils.Intersect(tetraA.Vertices[1], tetraA.Vertices[2], tetraA.Vertices[3], intersectCheck))
        {
            // flip23(A, B)
            Tetrahedron pabd = new(new[]
            {
                tetraA.Vertices[0], tetraA.Vertices[1], tetraA.Vertices[2], tetraB.Vertices[0],
            });
            Tetrahedron pbcd = new(new[]
            {
                tetraA.Vertices[0], tetraA.Vertices[2], tetraA.Vertices[3], tetraB.Vertices[0],
            });
            Tetrahedron pacd = new(new[]
            {
                tetraA.Vertices[0], tetraA.Vertices[1], tetraA.Vertices[3], tetraB.Vertices[0],
            });

            // push tetra pabd, pbcd, and pacd on stack
            newTetras.Enqueue(pabd);
            newTetras.Enqueue(pbcd);
            newTetras.Enqueue(pacd);
        }

        // case #2: concave && diagram(?) has tetra pdab
        else if (MathUtils.Intersect(tetraA.Vertices[1], tetraA.Vertices[2], tetraA.Vertices[3], intersectCheck)
                 && DelaunayDiagram.Contains(pdab))
        {
            // flip32(A, B, pdab)
            Tetrahedron pacd = new(new[]
            {
                tetraA.Vertices[0], tetraA.Vertices[1], tetraA.Vertices[3], tetraB.Vertices[0],
            });

            // push pacd and pdab on stack
            newTetras.Enqueue(pacd);
            newTetras.Enqueue(pdab);
        }

        float epsilon = MathUtils.EPSILON;
        for (int i = 0; i < 4; i++)
        {
            Triangle testTri = tetraB.Faces[i];
            bool coplanar = Orient(testTri.Vertices[0], testTri.Vertices[1], testTri.Vertices[2],
                    tetraA.Vertices[0])
                <= epsilon && Orient(testTri.Vertices[0], testTri.Vertices[1], testTri.Vertices[2],
                    tetraA.Vertices[0]) >= -epsilon;
            bool config44 = tetraA.Neighbors[i] == null;

            // case #3: degenerate coplanar && A and B are in config44 w/ C and D
            if (coplanar && config44)
            {
                // flip44(A, B, C, D)
                // push on stack the 4 tetra created
            }
        }

        /*
        // case #4: degenerate flat tetra
               flip23(A, B)
               push tetra pabd, pbcd, and pacd on stack
        */
    }

    private void DelIso(List<Tetrahedron> delaunay)
    {
        var voronoi = RecoverDiagram(delaunay);
        voronoi = RefineDiagram(voronoi);

        throw new NotImplementedException();
    }

    private ICollection<Cell> RefineDiagram(ICollection<Cell> diagram)
    {
        throw new NotImplementedException();
    }

    private ICollection<Cell> RecoverDiagram(List<Tetrahedron> delaunay)
    {
        throw new NotImplementedException();
    }

    /// <summary>
    ///   determines if a point p is over, under or lies on a plane defined by three points a, b and c
    /// </summary>
    /// <returns>
    ///   a positive value when the point p is above the plane defined by a, b and c; a negative value
    ///   if p is under the plane; and exactly 0 if p is directly on the plane.
    /// </returns>
    private float Orient(Vector3 a, Vector3 b, Vector3 c, Vector3 p)
    {
        float[,] orientMatrix =
        {
            { a.x, a.y, a.z, 1.0f },
            { b.x, b.y, b.z, 1.0f },
            { c.x, c.y, c.z, 1.0f },
            { p.x, p.y, p.z, 1.0f },
        };

        float determinant =
            orientMatrix[0, 0] * (orientMatrix[1, 1] * orientMatrix[2, 2] - orientMatrix[1, 2] * orientMatrix[2, 1]) -
            orientMatrix[0, 1] * (orientMatrix[1, 0] * orientMatrix[2, 2] - orientMatrix[1, 2] * orientMatrix[2, 0]) +
            orientMatrix[0, 2] * (orientMatrix[1, 0] * orientMatrix[2, 1] - orientMatrix[1, 1] * orientMatrix[2, 0]);

        return determinant;
    }

    /// <summary>
    ///   determines if a point p is inside, outside or lies on a sphere defined by four points a, b, c and d.
    /// </summary>
    /// <returns>
    ///   a positive value is returned if p is inside the sphere; a negative if p is outside; and exactly 0 if p
    ///   is directly on the sphere.
    /// </returns>
    private float InSphere(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3 p)
    {
        float[,] inSphereMatrix =
        {
            { a.x, a.y, a.z, (a.x * a.x) + (a.y * a.y) + (a.z * a.z), 1.0f },
            { b.x, b.y, b.z, (b.x * b.x) + (b.y * b.y) * (b.z * b.z), 1.0f },
            { c.x, c.y, c.z, (c.x * c.x) + (c.y * c.y) * (c.z * c.z), 1.0f },
            { d.x, d.y, d.z, (d.x * d.x) + (d.y * d.y) * (d.z * d.z), 1.0f },
            { p.x, p.y, p.z, (p.x * p.x) + (p.y * p.y) * (p.z * p.z), 1.0f },
        };

        float determinant =
            inSphereMatrix[0, 3] * inSphereMatrix[1, 2] * inSphereMatrix[2, 1] * inSphereMatrix[3, 0] -
            inSphereMatrix[0, 2] * inSphereMatrix[1, 3] * inSphereMatrix[2, 1] * inSphereMatrix[3, 0] -
            inSphereMatrix[0, 3] * inSphereMatrix[1, 1] * inSphereMatrix[2, 2] * inSphereMatrix[3, 0] +
            inSphereMatrix[0, 1] * inSphereMatrix[1, 3] * inSphereMatrix[2, 2] * inSphereMatrix[3, 0] +
            inSphereMatrix[0, 2] * inSphereMatrix[1, 1] * inSphereMatrix[2, 3] * inSphereMatrix[3, 0] -
            inSphereMatrix[0, 1] * inSphereMatrix[1, 2] * inSphereMatrix[2, 3] * inSphereMatrix[3, 0] -
            inSphereMatrix[0, 3] * inSphereMatrix[1, 2] * inSphereMatrix[2, 0] * inSphereMatrix[3, 1] +
            inSphereMatrix[0, 2] * inSphereMatrix[1, 3] * inSphereMatrix[2, 0] * inSphereMatrix[3, 1] +
            inSphereMatrix[0, 3] * inSphereMatrix[1, 0] * inSphereMatrix[2, 2] * inSphereMatrix[3, 1] -
            inSphereMatrix[0, 0] * inSphereMatrix[1, 3] * inSphereMatrix[2, 2] * inSphereMatrix[3, 1] -
            inSphereMatrix[0, 2] * inSphereMatrix[1, 0] * inSphereMatrix[2, 3] * inSphereMatrix[3, 1] +
            inSphereMatrix[0, 0] * inSphereMatrix[1, 2] * inSphereMatrix[2, 3] * inSphereMatrix[3, 1] +
            inSphereMatrix[0, 3] * inSphereMatrix[1, 1] * inSphereMatrix[2, 0] * inSphereMatrix[3, 2] -
            inSphereMatrix[0, 1] * inSphereMatrix[1, 3] * inSphereMatrix[2, 0] * inSphereMatrix[3, 2] -
            inSphereMatrix[0, 3] * inSphereMatrix[1, 0] * inSphereMatrix[2, 1] * inSphereMatrix[3, 2] +
            inSphereMatrix[0, 0] * inSphereMatrix[1, 3] * inSphereMatrix[2, 1] * inSphereMatrix[3, 2] +
            inSphereMatrix[0, 1] * inSphereMatrix[1, 0] * inSphereMatrix[2, 3] * inSphereMatrix[3, 2] -
            inSphereMatrix[0, 0] * inSphereMatrix[1, 1] * inSphereMatrix[2, 3] * inSphereMatrix[3, 2] -
            inSphereMatrix[0, 2] * inSphereMatrix[1, 1] * inSphereMatrix[2, 0] * inSphereMatrix[3, 3] +
            inSphereMatrix[0, 1] * inSphereMatrix[1, 2] * inSphereMatrix[2, 0] * inSphereMatrix[3, 3] +
            inSphereMatrix[0, 2] * inSphereMatrix[1, 0] * inSphereMatrix[2, 1] * inSphereMatrix[3, 3] -
            inSphereMatrix[0, 0] * inSphereMatrix[1, 2] * inSphereMatrix[2, 1] * inSphereMatrix[3, 3] -
            inSphereMatrix[0, 1] * inSphereMatrix[1, 0] * inSphereMatrix[2, 2] * inSphereMatrix[3, 3] +
            inSphereMatrix[0, 0] * inSphereMatrix[1, 1] * inSphereMatrix[2, 2] * inSphereMatrix[3, 3];

        return determinant;
    }

    // voronoi cell
    public struct Cell
    {
        public List<Vector3> Vertices;
        public Vector3 Site;

        public Cell(Vector3 seed, List<Vector3> verts)
        {
            Site = seed;
            Vertices = verts;
        }
    }

    // delaunay triangle
    public struct Triangle
    {
        public Vector3[] Vertices;

        public Triangle(float[] a, float[] b, float[] c)
        {
            Vertices = new[]
            {
                new Vector3(a[0], a[1], a[2]),
                new Vector3(b[0], b[1], b[2]),
                new Vector3(c[0], c[1], c[2]),
            };
        }
    }

    // a tetrahedron that helps form a 3D delaunay diagram
    public struct Tetrahedron : IEquatable<Tetrahedron>
    {
        public Vector3[] Vertices;
        public Triangle[] Faces;
        public Tetrahedron*[] Neighbors;

        public Tetrahedron(Vector3[] verts)
        {
            Vertices = verts;

            float[] q = { verts[0].x, verts[0].y, verts[0].z };
            float[] r = { verts[1].x, verts[1].y, verts[1].z };
            float[] l = { verts[2].x, verts[2].y, verts[2].z };
            float[] t = { verts[3].x, verts[3].y, verts[3].z };

            Faces = new[]
            {
                new Triangle(q, r, l),
                new Triangle(q, t, r),
                new Triangle(q, l, t),
                new Triangle(r, t, l),
            };

            // edge flags to ID neighbors; each int is the index of the neighbor
            Neighbors = new Tetrahedron*[4];
        }

        public static bool operator ==(Tetrahedron tetra, Tetrahedron other)
        {
            return tetra.Equals(other);
        }

        public static bool operator !=(Tetrahedron tetra, Tetrahedron other)
        {
            return !tetra.Equals(other);
        }

        public bool Equals(Tetrahedron other)
        {
            return Vertices == other.Vertices;
        }

        public override bool Equals(object obj)
        {
            if (!(obj is Tetrahedron tetra))
                return false;

            return Equals(tetra);
        }

        public override int GetHashCode()
        {
            var hashCode = -1949845991;

            hashCode = hashCode * -1287342897 + Vertices.GetHashCode();
            hashCode = hashCode * -1287342897 + Faces.GetHashCode();
            hashCode = hashCode * -1287342897 + Neighbors.GetHashCode();

            return hashCode;
        }
    }
}
