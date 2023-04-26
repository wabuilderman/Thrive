using System;
using System.Collections.Generic;
using System.Numerics;

public class VoxelGrid
{
    public int Width;
    public int Length;
    public int Height;

    public List<int> Voxels;

    private int isoLevel = 1;

    public VoxelGrid(int x, int y, int z)
    {
        Width = x;
        Height = y;
        Length = z;

        Voxels = new List<int>(x * y * z);
    }

    public ICollection<Vector3> GetRandomSeeds(MetaballLayout<Metaball> fields)
    {
        Vector3 fieldA = default(Vector3);
        Vector3 fieldB = default(Vector3);

        Random rng = new();
        var activeVoxels = new List<Vector3>();

        for (int i = 0; i < fields.Count; i++)
        {
            fieldA.X = fields[i].Position.x;
            fieldA.Y = fields[i].Position.y;
            fieldA.Z = fields[i].Position.z;

            fieldB.X = fields[i + 1].Position.x;
            fieldB.Y = fields[i + 1].Position.y;
            fieldB.Z = fields[i + 1].Position.z;

            int point = Voxels.Random(rng);
            var voxelPos = new Vector3(point / (Width * Height),
                (point - Length * Width * Height) / Width,
                point % Width);

            bool voxelActive = Math.Abs(ConvolutionSurface.GetIntegralAtPoint(fieldA, fieldB, voxelPos, 1.0f)
                - isoLevel) < 1.0f;

            if (voxelActive)
                activeVoxels.Add(voxelPos);
        }

        ICollection<Vector3> seeds = activeVoxels;

        return seeds;
    }
}
