using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using Vector3 = Godot.Vector3;

public class Voronoi
{
    /*
     * seed - points where cells are generated
     * cell - domain of all points closer to site[n] than any other site
     * del - Delaunay dual of the voronoi diagram
     */

    public List<Cell> Diagram;
    public Voronoi(List<Vector3> seeds)
    {
        Diagram = GetCellBounds(seeds);
    }

    private List<Cell> GetCellBounds(List<Vector3> sites)
    {
        List<Cell> voronoiCells = null;
        var slice = new List<Vector3>();

        var sliceQueue = new ConcurrentQueue<List<Vector3>>();
        for (int i = 0; i < (int)sites[sites.Count - 1].z; i++)
        {
            foreach (var site in sites)
            {
                if ((int)site.z == i)
                    slice.Add(site);
            }

            sliceQueue.Enqueue(slice);
        }

        Action map = () =>
        {
            // point between cells with a line dividing them
            List<Vector3> relations;
            while (sliceQueue.TryDequeue(out relations))
            {
            }
        };

        // try to map multiple slices at once
        Parallel.Invoke(map, map, map, map);
        return voronoiCells;
    }

    public struct Cell
    {
        // public List<Vector3> Corners;
        public Vector3 Site;

        public Cell(Vector3 seed)
        {
            Site = seed;
        }
    }
}
