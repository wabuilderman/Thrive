using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using Godot.Collections;
using Supercluster.KDTree;
using Array = Godot.Collections.Array;
using Vector3 = Godot.Vector3;

public class Voronoi
{
    public ICollection<Cell>? VoronoiDiagram;
    public ICollection<Simplex>? DelaunayDiagram;

    // very big tetrahedron, will automate making this
    private readonly float[][] bigTetra =
    {
        new[] { 0, -24999.75f, -74999.25f },
        new[] { 0, 74999.25f, 0 },
        new[] { 53032.61968f, -24999.75f, -53032.61968f },
        new[] { -53032.61968f, -24999.75f, -53032.61968f },
    };

    // standard euclidean distance
    private readonly Func<float[], float[], double> l2Norm = (p, q) =>
    {
        double dist = 0;
        for (int i = 0; i < p.Length; i++)
        {
            dist += (p[i] - q[i]) * (p[i] - q[i]);
        }

        return dist;
    };

    public Voronoi(List<Vector3> seeds)
    {
        DelaunayDiagram = InitializeDiagram(seeds);
    }

    /// <summary>
    /// determines if a point p is over, under or lies on a plane defined by three points a, b and c
    /// </summary>
    /// <returns>
    /// a positive value when the point p is above the plane defined by a, b and c; a negative value
    /// if p is under the plane; and exactly 0 if p is directly on the plane.
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
    /// determines if a point p is inside, outside or lies on a sphere defined by four points a, b, c and d.
    /// </summary>
    /// <returns>
    /// a positive value is returned if p is inside the sphere; a negative if p is outside; and exactly 0 if p
    /// is directly on the sphere.
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

    /// <summary>
    ///     <para>
    ///         O.Devillers, S.Pion, and M.Teillaud 'Walking in a Triangulation'
    ///     </para>
    /// International Journal of Foundations of Computer Science, 13(2):181-199, 2002
    ///     <para>
    ///         Available from: https://inria.hal.science/inria-00102194/document
    ///     </para>
    /// </summary>
    /// <param name="tet">tetrahedron to walk from</param>
    /// <param name="point">query point</param>
    private Tetrahedron Walk(Tetrahedron tet, Vector3 point)
    {
        // Remembering Stochastic Walk
        // from q(uery) to p(oint). tet=abcd is a tetrahedron of simplexes (faces)
        // previous = tet; end = false;
        var previous = tet;
        bool end = false;

        while (!end)
        {
            // face = random facet of tet;
            int random = new Random().Next(3);
            Simplex face = tet.Faces[random];
            bool neighbor = false;
            bool otherSide = false;

            // if( point not neighbor of previous through face ) && ( point on other side of face )
            if (!neighbor && otherSide)
            {
                // previous = tet;
                // tet = neighbor( tet through face);
            }
            else
            {
                // face = next facet of tet;
                face = random < 4 ? tet.Faces[random++] : tet.Faces[0];

                // if( point not neighbor of previous through face ) && ( point on other side of face )
                if (!neighbor && otherSide)
                {
                    // previous = tet;
                    // tet = neighbor( tet through face);
                }
                else
                {
                    end = true;
                }
            }
        }

        // t contains p

        throw new NotImplementedException();
    }

    private void Flip()
    {
        throw new NotImplementedException();
    }

    // this inserts a point and calculates new tetrahedrons
    private void InsertPoint(ICollection<Simplex> delaunay, Vector3 point)
    {
        // tet <--walk
        // Tetrahedron tet = Walk(delaunay, point);

        // insert point in tet with a flip14
        // push 4 new tets on stack
        // while stack is non-empty do
        //   tet = {p, a, b, c} <--pop from stack
        //   tet[a] = {a, b, c, d} <--get adjacent tet of delaunay having abc as facet
        //   if d is inside circumsphere of delaunay then
        //       Flip(delaunay, delaunay[a])
        //   end if
        // end while

        throw new NotImplementedException();
    }

    /// <summary>
    /// <para>Hugo Ledoux 'Computing the 3D Voronoi Diagram Robustly: An Easy Explanation'</para>
    /// Delft University of Technology (OTB-section GIS Technology) [Internet] 2007
    /// <para>Available from: http://www.gdmc.nl/publications/2007/Computing_3D_Voronoi_Diagram.pdf</para>
    /// </summary>
    /// <returns>
    /// initialized voronoi diagram
    /// </returns>
    private List<Simplex>? InitializeDiagram(List<Vector3> seeds)
    {
        // big tetrahedron that contains all points to start
        var triangulation = new List<Simplex>();
        triangulation.Add(new Simplex(bigTetra[0], bigTetra[2], bigTetra[1]));
        triangulation.Add(new Simplex(bigTetra[0], bigTetra[3], bigTetra[1]));
        triangulation.Add(new Simplex(bigTetra[0], bigTetra[2], bigTetra[3]));
        triangulation.Add(new Simplex(bigTetra[1], bigTetra[2], bigTetra[3]));

        // insert seeds as query points, then rebuild diagram until we triangulate every point.
        // gives low poly initial triangulation that we make more detailed with Del-Iso
        for (int i = 0; i < seeds.Count; i++)
        {
            Vector3 query = seeds[i];
        }

        return triangulation;
    }

    private ICollection<Cell> RefineDiagram(ICollection<Cell> diagram)
    {
        throw new NotImplementedException();
    }

    private ICollection<Cell> RecoverDiagram(ICollection<Cell> diagram)
    {
        throw new NotImplementedException();
    }

    // find nearest neighbours using a kd tree
    private IEnumerable<Vector3> GetNearestNeighbours(List<Vector3> sites)
    {
        int siteCount = sites.Count;

        var data = new List<float[]>();

        for (int i = 0; i < siteCount; i++)
        {
            data.Add(new float[] { sites[i].x, sites[i].y, sites[i].z });
        }

        float[][] treeData = data.ToArray();
        var treeNodes = sites.Select(p => p.ToString()).ToArray();
        var tree = new KDTree<float, string>(3, treeData, treeNodes, l2Norm);

        // gotta work on this one a lot
        throw new NotImplementedException();
    }

    // voronoi cell
    public struct Cell
    {
        public List<Vector3> Corners;
        public Vector3 Site;

        public Cell(Vector3 seed, List<Vector3> verts)
        {
            Site = seed;
            Corners = verts;
        }
    }

    // delaunay triangle
    public struct Simplex
    {
        public Array<Vector3> Corners;
        public Vector3 Centroid;

        public Simplex(float[] a, float[] b, float[] c)
        {
            Corners = new Array<Vector3>
            {
                new(a[0], a[1], a[2]),
                new(b[0], b[1], b[2]),
                new(c[0], c[1], c[2]),
            };

            float centroidX = (Corners[0].x + Corners[1].x + Corners[2].x) / 3;
            float centroidY = (Corners[0].y + Corners[1].y + Corners[2].y) / 3;
            float centroidZ = (Corners[0].z + Corners[1].z + Corners[2].z) / 3;
            Centroid = new Vector3(centroidX, centroidY, centroidZ);
        }
    }

    // a tetrahedron that helps form a 3D delaunay diagram
    public struct Tetrahedron
    {
        public Array<Simplex> Faces;

        public Tetrahedron(Array<Simplex> facets)
        {
            Faces = facets;
        }
    }
}
