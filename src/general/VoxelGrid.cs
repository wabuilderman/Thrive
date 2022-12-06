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

    public void GetActiveVoxels(MetaballLayout<Metaball> fields)
    {
        Random rng = new Random();
        int point = Voxels.Random(rng);
        Vector3 voxelPos = new Vector3(point / (Width * Height),
            (point - Length * Width * Height) / Width,
            point % Width);

        bool voxelActive;
        Vector3 fieldA = default(Vector3);
        Vector3 fieldB = default(Vector3);
        for (int i = 0; i < fields.Count; i++)
        {
            fieldA.X = fields[i].Position.x;
            fieldA.Y = fields[i].Position.y;
            fieldA.Z = fields[i].Position.z;

            fieldB.X = fields[i + 1].Position.x;
            fieldB.Y = fields[i + 1].Position.y;
            fieldB.Z = fields[i + 1].Position.z;

            voxelActive = Math.Abs(ConvolutionSurface.GetIntegralAtPoint(fieldA, fieldB, voxelPos, 1.0f)
                - isoLevel) < MathUtils.EPSILON;
        }
    }
}
