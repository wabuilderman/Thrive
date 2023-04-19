using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using Supercluster.KDTree;
using Vector3 = Godot.Vector3;

public class Voronoi
{
    public ICollection<Cell>? Diagram;

    // standard euclidean distance
    private Func<float[], float[], double> l2Norm = (p, q) =>
    {
        double dist = 0;
        for (int i = 0; i < p.Length; i++)
        {
            dist += (p[i] - q[i]) * (p[i] - q[i]);
        }

        return dist;
    };

    // very big tetrahedron, will automate making this
    private float[,] bigTetra =
    {
        { 0, -24999.75f, -74999.25f },
        { 0, 74999.25f, 0 },
        { 53032.61968f, -74999.25f, 53032.61968f },
        { -53032.61968f, -74999.25f, -53032.61968f },
    };

    public Voronoi(List<Vector3> seeds)
    {
        Diagram = InitializeDiagram(seeds);
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
    /// <para>Hugo Ledoux 'Computing the 3D Voronoi Diagram Robustly: An Easy Explanation'</para>
    /// Delft University of Technology (OTB-section GIS Technology) [Internet] 2007
    /// <para>Available from: http://www.gdmc.nl/publications/2007/Computing_3D_Voronoi_Diagram.pdf</para>
    /// </summary>
    /// <returns>
    /// initialized voronoi diagram
    /// </returns>
    private List<Cell>? InitializeDiagram(List<Vector3> seeds)
    {
        var voronoiCells = new List<Cell>();

        for (int i = 0; i < seeds.Count; i++)
        {
            Vector3 root = seeds[i];
        }

        return voronoiCells;
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
        return null;
    }

    public struct Cell
    {
        public ICollection<Vector3> Corners;
        public Vector3 Site;

        public Cell(Vector3 seed, List<Vector3> verts)
        {
            Site = seed;
            Corners = verts;
        }
    }
}
